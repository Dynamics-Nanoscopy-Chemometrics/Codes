% ========================================================================
% AslsBaseSingle Classic ASLS baseline for a single vector (ASLS)
% ========================================================================
% This code implements the asymmetric least squares (ASLS) baseline
% correction method originally proposed by:
%
%   Eilers, P. H. C., & Boelens, H. F. M. (2005).
%   Baseline correction with asymmetric least squares smoothing.
%   Leiden University Medical Centre Report, 1(1), 5.
%
% ------------------------------
% HOW TO CITE (per Lovelace's Square)
% ------------------------------
% If you use this code, please cite BOTH:
%
% (1) The original method:
%     Eilers, P. H. C., & Boelens, H. F. M. (2005).
%     Baseline correction with asymmetric least squares smoothing.
%     Leiden University Medical Centre Report, 1(1), 5.
%
% (2) This implementation (example format; update the URL when available):
%     Gómez, Y. R. (2025). AslsBaseSingle ASLS baseline correction (MATLAB).
%     Lovelace's Square. Version 1.0.0. Last accessed: 2025-11-12.
%
% ------------------------------
% REFERENCES (BibTeX)
% ------------------------------
% @article{Eilers2005AsLS,
%   author  = {Eilers, Paul H. C. and Boelens, Hans F. M.},
%   title   = {Baseline Correction with Asymmetric Least Squares Smoothing},
%   journal = {Leiden University Medical Centre Report},
%   year    = {2005},
%   volume  = {1},
%   number  = {1},
%   pages   = {5}
% }
%
% ------------------------------
% METADATA
% ------------------------------
% Authors: Yesid Roman Gomez
% Date Created: 2025-10-16
% License: MIT
% Version: 1.0.0
% Implements: Eilers & Boelens (2005) - non-novel implementation
% Reviewed by Lovelace's Square team: No
%
% ------------------------------
% USAGE
% ------------------------------
% Args:
%   y (nx1 double): signal (column or row vector accepted)
%   lambda (double): smoothness penalty (>= 0)
%   p (double): asymmetry parameter in [0, 1]
%
% Returns:
%   z (nx1 double): estimated baseline
%   Wout (nxk double): weight history (k iterations)
%
% Example:
%   z = AslsBaseSingle(y, 1e5, 0.001);
% ========================================================================

function [z, Wout] = AslsBaseSingle(y, lambda, p) 
% ASLSBASESINGLE Classic ASLS baseline for a single vector

    y = y(:);
    m = numel(y);

    if m == 0
        z = y;
        Wout = [];
        return;
    end

    % Parameters / constants
    MAX_ITERS = 10;

    % Precompute 2nd-difference operator
    D  = diff(speye(m), 2);
    w  = ones(m, 1);
    z  = zeros(m, 1);
    Wout = [];

    Ieps = 1e-12 * speye(m);
    DtD  = D' * D;

    for i = 1:MAX_ITERS
        W = spdiags(w, 0, m, m);
        % (W + lambda * D' * D) z = W * y  via Cholesky
        C = chol(W + lambda * DtD + Ieps, 'lower');
        z = C' \ (C \ (w .* y));

        % Update asymmetric weights
        w = p * (y > z) + (1 - p) * (y < z);

        % Track weights (optional)
        Wout = [Wout, w]; %#ok<AGROW>
    end
end

