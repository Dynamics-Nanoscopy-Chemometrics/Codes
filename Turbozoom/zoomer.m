function [x, it, err_hist] = zoomer(img, k, lambda, delta, zoomf, nit, sigma)
%ZOOMER  Single-image super-resolution by quadratic regularized CG.
%
%   Reference:
%   Eilers, P. H. C., & Ruckebusch, C. (2022). Fast and simple
%   super-resolution with single images. Scientific Reports, 12, 11241.
%   https://doi.org/10.1038/s41598-022-14874-8
%
%   Paul H. C. Eilers - Dept. of Biostatistics, Erasmus MC, Rotterdam, NL
%   Cyril Ruckebusch - LASIRE CNRS, University of Lille, Lille, France
%
%   Version: 1.1
%   Date:    November 2025
%   License: MIT
%
%   [X, IT, ERR_HIST] = ZOOMER(IMG, K, LAMBDA, DELTA, ZOOMF, NIT, SIGMA)
%   performs image super-resolution by solving a regularized least-squares
%   problem with a conjugate gradient (CG) scheme.
%
%   The forward model is
%
%       IMG ≈ S1 * X * S2'
%
%   where:
%       - X   is the unknown high-resolution image,
%       - IMG is the observed low-resolution image,
%       - S1, S2 combine high-resolution Gaussian blur and binning
%         (downsampling) along each dimension.
%
%   The optimization problem is:
%
%       min_X  ||S1*X*S2' - IMG||_F^2                                   ...
%              + k * ||X||_F^2                                         ...
%              + lambda * ( ||D1*X||_F^2 + ||X*D2'||_F^2 )
%
%   where:
%       - k      is a zeroth-order Tikhonov (ridge) parameter on X,
%       - lambda controls first-order Tikhonov (smoothness) penalties
%         via discrete difference operators D1 and D2.
%
% Inputs:
%   img     - Low-resolution input image (n1 x n2 matrix) to be
%             super-resolved.
%   k       - Zeroth-order Tikhonov (ridge) regularization parameter on X
%             (scalar, k >= 0).
%   lambda  - First-order Tikhonov regularization parameter controlling
%             smoothness through finite differences in rows and columns
%             (scalar, lambda >= 0).
%   delta   - Convergence threshold on the CG residual RMS
%             (scalar, delta > 0).
%   zoomf   - Zoom factor (p1 = zoomf*n1, p2 = zoomf*n2); must be a
%             positive integer (e.g., 2 for 2× zoom).
%   nit     - Maximum number of CG iterations (positive integer).
%   sigma   - Standard deviation for the Gaussian PSF used to construct
%             the blur operators (scalar, sigma > 0).
%
% Outputs:
%   x        - Super-resolved image (p1 x p2 matrix, with p1 = zoomf*n1
%              and p2 = zoomf*n2).
%   it       - Number of iterations performed before convergence or
%              reaching nit.
%   err_hist - Column vector of length IT containing the data-space
%              RMS error at each iteration:
%
%                  err_hist(it) = ||S1*X_it*S2' - IMG||_F / sqrt(numel(IMG))
%
% Algorithm (high level):
%   1. Construct 1D Gaussian PSF matrices and corresponding blur+bin
%      operators S1, S2.
%   2. Form the normal-equation components:
%          G1 = S1'*S1,   G2 = S2'*S2,   U = S1'*IMG*S2.
%   3. Solve the quadratic problem with CG on X:
%          A(X) = G1*X*G2 + k*X + V1*X + X*V2',
%      where V1, V2 encode the smoothness (difference) penalties.
%   4. Track both the CG residual RMS (for stopping) and the data-space
%      RMS error via the full forward model S1*X*S2' - IMG.
%
% Example:
%   img = imread('lowres.png');
%   img = rgb2gray(img);
%   zoomf = 2;
%   [x, it, err_hist] = zoomer(img, 1e-2, 1e-1, 1e-6, zoomf, 100, 1.5);
%   imshow(x, []);
%
% See also: IMRESIZE, DECONVLUCY, PSFEST
%% Input validation
if nargin < 7
    error('zoomer:notEnoughInputs', 'All 7 input arguments are required.');
end

if zoomf <= 0 || floor(zoomf) ~= zoomf
    error('zoomer:invalidZoomFactor', 'Zoom factor must be a positive integer.');
end

if nit <= 0 || floor(nit) ~= nit
    error('zoomer:invalidIterations', 'Number of iterations must be a positive integer.');
end

if lambda < 0 || k < 0 || delta <= 0 || sigma <= 0
    error('zoomer:invalidParameters', ...
        'Parameters lambda, k, sigma must be non-negative and delta and sigma must be positive.');
end

%% PSF (Point Spread Function) Estimation and blur+bin operators
[n1, n2] = size(img);      % Low-resolution image dimensions
p1 = zoomf * n1;           % High-resolution height
p2 = zoomf * n2;           % High-resolution width

% 1D PSF matrices (high-resolution)
s1a = psfest(p1, sigma);   % PSF for dimension 1 (p1 x p1)
s2a = psfest(p2, sigma);   % PSF for dimension 2 (p2 x p2)

% Blur + binning operators:
% kron(eye(n1), ones(1, zoomf)) implements block summation (binning).
s1 = kron(eye(n1), ones(1, zoomf)) * s1a;  % (n1 x p1)
s2 = kron(eye(n2), ones(1, zoomf)) * s2a;  % (n2 x p2)

%% Normal-equation components (data term)
% Forward model: IMG ≈ S1 * X * S2'
% Normal equations: G1*X*G2 + ... = U
u  = s1' * img * s2;       % Right-hand side U (p1 x p2)
g1 = s1' * s1;             % G1 (p1 x p1)
g2 = s2' * s2;             % G2 (p2 x p2)

%% Conjugate Gradient Setup
n1_zoom = size(g1, 2);     % = p1
n2_zoom = size(g2, 2);     % = p2
x = zeros(n1_zoom, n2_zoom);   % Initial solution (high-resolution)

% Initial residual r0 = b - A*x0. With x0 = 0, r0 = U.
r = u;
p = r;

% Preallocate data-space error history
err_hist = zeros(nit, 1);

% Difference operators for smoothness (first-order Tikhonov)
d1 = diff(eye(n1_zoom));   % (p1-1 x p1)
d2 = diff(eye(n2_zoom));   % (p2-1 x p2)

% Smoothness penalty matrices
v1 = lambda * (d1' * d1);  % V1 (p1 x p1)
v2 = lambda * (d2' * d2);  % V2 (p2 x p2)

%% Conjugate Gradient Iteration
for it = 1:nit
    % Apply the system matrix A to p:
    % A(p) = G1*p*G2 + k*p + V1*p + p*V2'
    q = g1 * p * g2 + k * p + v1 * p + p * v2';

    % Step size alpha = (r,r) / (p, A p)
    rs1 = sum(r(:).^2);        % ||r||^2
    denom = sum(p(:) .* q(:));
    if denom == 0
        % Degenerate direction; stop early.
        err_hist = err_hist(1:it-1);
        break;
    end
    alpha = rs1 / denom;

    % Update solution
    x = x + alpha * p;
    
    % Update residual
    rn = r - alpha * q;
    rs2 = sum(rn(:).^2);       % ||r_new||^2

    % CG coefficient beta
    beta = rs2 / rs1;

    % Update search direction and residual
    p = rn + beta * p;
    r = rn;

    % CG residual RMS (used for stopping criterion)
    rms_cg = sqrt(rs2 / (n1_zoom * n2_zoom));

    % Data-space error (includes blur + binning: S1 * X * S2' - IMG)
    diff_img      = s1 * x * s2' - img;           % low-res domain
    err_hist(it)  = sqrt(mean(diff_img(:).^2));   % RMS in measurement space

    % Check convergence
    if rms_cg < delta
        err_hist = err_hist(1:it);   % Trim unused entries
        break;
    end
end

end
