% TEST  Demonstration script for the ZOOMER super-resolution algorithm.
%
%   This script:
%       1) Generates a synthetic high-resolution image made of:
%            - finite-size emitters arranged in a circle or figure-8,
%            - optional filament-like continuous lines with controllable
%              length and thickness.
%       2) Simulates a low-resolution measurement by applying Gaussian blur
%          and binning (downsampling).
%       3) Uses the ZOOMER function to reconstruct a high-resolution image
%          from the low-resolution data.
%       4) Displays three images side by side:
%              - Ground truth (high-resolution)
%              - Low-resolution blurred image
%              - Reconstructed high-resolution image
%
%   Requirements on the MATLAB path:
%       - zoomer.m   : Super-resolution reconstruction function
%       - psfest.m   : 1D PSF matrix generator used by zoomer
%       - binning.m  : Binning (downsampling) helper function
%
%   Date created: 2025-11-12
%   License: MIT
%   Version: 1.1  (added filament-like continuous lines)
%
%   Reference for the method:
%       Eilers, P. H. C., & Ruckebusch, C. (2022).
%       Fast and simple super-resolution with single images.
%       Scientific Reports, 12, 11241.

%% Clean workspace and set reproducible random seed
clearvars;
close all;
clc;

rng(1, "twister");   % Reproducible randomness for emitter positions/intensities

fprintf("=== Zoomer Super-Resolution Demo ===\n\n");

%% 1. Configuration and synthetic high-resolution image
fprintf("Step 1: Creating synthetic high-resolution image...\n");

% ----- Global configuration (treated as constants) -----------------------
HR_SIZE      = 200;        % High-resolution image size (square)
ZOOM_FACTOR  = 3;          % Downsampling / reconstruction zoom factor

% Discrete emitter configuration (circle or figure-8)
SHAPE_TYPE   = "figure8";  % "circle" or "figure8"
NUM_EMITTERS = 50;         % Number of emitters along the shape

EMITTER_SHAPE  = "circle"; % "pixel", "square3", "square4", "circle"
EMITTER_RADIUS = 3;        % Radius (for circle) or half-size for squares

% Filament configuration (continuous line structures)
ADD_FILAMENTS        = true;   % Set false if you want to disable filaments
NUM_FILAMENTS        = 6;      % Number of filaments
FILAMENT_LENGTH_PIX  = 200;    % Approximate length of each filament (in pixels)
FILAMENT_THICKNESS   = 3;      % Full thickness (diameter) in pixels

% Forward-simulation parameters
SIGMA_PSF_SIM  = 3.0;   % PSF sigma used for forward simulation
NOISE_STD      = 0.0025; % Set > 0 to add Gaussian noise to LR image

% Zoom and image sizes
highResSize = HR_SIZE;
zoomFactor  = ZOOM_FACTOR;
lowResSize  = highResSize / zoomFactor;

% Preallocate high-resolution image
imageHighRes = zeros(highResSize, highResSize);

% Precompute coordinate grid (used for circles and filaments)
[xxGrid, yyGrid] = meshgrid(1:highResSize, 1:highResSize);

% ----- Generate emitter positions along chosen shape --------------------
centerX = highResSize / 2;
centerY = highResSize / 2;

switch lower(SHAPE_TYPE)
    case "circle"
        radius = 70;
        theta  = linspace(0, 2 * pi, NUM_EMITTERS + 1);
        theta(end) = [];

        xPositions = centerX + radius * cos(theta);
        yPositions = centerY + radius * sin(theta);

    case "figure8"
        radius = 40;
        offset = 50;

        nEach = max(1, floor(NUM_EMITTERS / 2));
        theta = linspace(0, 2 * pi, nEach + 1);
        theta(end) = [];

        xLeft  = centerX - offset / 2 + radius * cos(theta);
        yLeft  = centerY + radius * sin(theta);

        xRight = centerX + offset / 2 + radius * cos(theta);
        yRight = centerY - radius * sin(theta);

        xPositions = [xLeft,  xRight];
        yPositions = [yLeft,  yRight];

    otherwise
        error("Unknown SHAPE_TYPE: %s (use ""circle"" or ""figure8"").", ...
            SHAPE_TYPE);
end

% Add small jitter to positions to simulate imperfect placement
positionJitter = 0.5;
xPositions = xPositions + randn(size(xPositions)) * positionJitter;
yPositions = yPositions + randn(size(yPositions)) * positionJitter;

% Round and clamp to valid pixel indices
xPositions = round(min(max(xPositions, 1), highResSize));
yPositions = round(min(max(yPositions, 1), highResSize));

% ----- Paint finite-size emitters on the high-resolution grid ----------
for emitterIndex = 1:numel(xPositions)
    x = xPositions(emitterIndex);
    y = yPositions(emitterIndex);

    % Each emitter has slightly different intensity
    intensity = 0.8 + 0.4 * rand();

    switch lower(EMITTER_SHAPE)
        case "pixel"
            imageHighRes(y, x) = imageHighRes(y, x) + intensity;

        case "square3"
            halfSize = 1;  % 3 x 3
            yRange = max(1, y - halfSize) : min(highResSize, y + halfSize);
            xRange = max(1, x - halfSize) : min(highResSize, x + halfSize);
            imageHighRes(yRange, xRange) = ...
                imageHighRes(yRange, xRange) + intensity;

        case "square4"
            halfSize = 2;  % 4 x 4
            yRange = max(1, y - halfSize + 1) : ...
                     min(highResSize, y + halfSize);
            xRange = max(1, x - halfSize + 1) : ...
                     min(highResSize, x + halfSize);
            imageHighRes(yRange, xRange) = ...
                imageHighRes(yRange, xRange) + intensity;

        case "circle"
            mask = ((xxGrid - x).^2 + (yyGrid - y).^2) <= EMITTER_RADIUS^2;
            imageHighRes(mask) = imageHighRes(mask) + intensity;

        otherwise
            error("Unknown EMITTER_SHAPE: %s", EMITTER_SHAPE);
    end
end

fprintf("   - Created %d emitters (%s, %s)\n", ...
    numel(xPositions), SHAPE_TYPE, EMITTER_SHAPE);

% ----- Add filament-like continuous lines (optional) --------------------
if ADD_FILAMENTS
    fprintf("   - Adding %d filament(s) (length ≈ %d px, thickness ≈ %d px)\n", ...
        NUM_FILAMENTS, FILAMENT_LENGTH_PIX, FILAMENT_THICKNESS);

    halfThickness = FILAMENT_THICKNESS / 2;
    radiusSqFil   = (halfThickness)^2;

    for fIdx = 1:NUM_FILAMENTS
        % Random center near the middle of the image
        cx = centerX + 10 * randn();
        cy = centerY + 10 * randn();

        % Random orientation
        angle = 2 * pi * rand();
        halfLen = FILAMENT_LENGTH_PIX / 2;

        % Endpoints of the filament before clamping
        x0 = cx - halfLen * cos(angle);
        y0 = cy - halfLen * sin(angle);
        x1 = cx + halfLen * cos(angle);
        y1 = cy + halfLen * sin(angle);

        % Clamp endpoints to image bounds
        x0 = min(max(x0, 1), highResSize);
        y0 = min(max(y0, 1), highResSize);
        x1 = min(max(x1, 1), highResSize);
        y1 = min(max(y1, 1), highResSize);

        % Direction vector of the segment
        vx = x1 - x0;
        vy = y1 - y0;
        denom = vx^2 + vy^2;
        if denom < eps
            % Degenerate filament (too short after clamping); skip it.
            continue;
        end

        % Compute projection of each pixel onto the segment [x0,y0]-[x1,y1]
        % Parameter t in [0,1] along the segment
        t = ((xxGrid - x0) * vx + (yyGrid - y0) * vy) / denom;
        t = max(0, min(1, t));  % clamp to the segment

        % Closest point on the segment for each pixel
        xClosest = x0 + t * vx;
        yClosest = y0 + t * vy;

        % Squared distance from each pixel to the filament
        distSq = (xxGrid - xClosest).^2 + (yyGrid - yClosest).^2;

        % Mask for pixels within the filament thickness
        filamentMask = distSq <= radiusSqFil;

        % Filament intensity (can be tuned; here slightly lower than emitters)
        filamentIntensity = 0.7 + 0.3 * rand();

        % Add filament contribution
        imageHighRes(filamentMask) = ...
            imageHighRes(filamentMask) + filamentIntensity;
    end
else
    fprintf("   - Filaments disabled (ADD_FILAMENTS = false)\n");
end

fprintf("   - High-resolution size: %d x %d pixels\n", ...
    highResSize, highResSize);

%% 2. Simulate low-resolution blurred image (PSF + binning + noise)
fprintf("\nStep 2: Simulating low-resolution measurement...\n");

% ----- Build Gaussian PSF for forward simulation -----------------------
psfHalfWidth = ceil(3 * SIGMA_PSF_SIM);
psfSize      = 2 * psfHalfWidth + 1;

[xGridPSF, yGridPSF] = meshgrid(-psfHalfWidth:psfHalfWidth, ...
                                -psfHalfWidth:psfHalfWidth);

psfKernel = exp(-(xGridPSF.^2 + yGridPSF.^2) ./ ...
                (2 * SIGMA_PSF_SIM^2));

psfKernel = psfKernel / sum(psfKernel(:));

% Convolve high-resolution image with PSF
imageHighResBlurred = conv2(imageHighRes, psfKernel, "same");

fprintf("   - Applied Gaussian PSF (sigma = %.2f)\n", SIGMA_PSF_SIM);
fprintf("   - PSF kernel size: %d x %d\n", psfSize, psfSize);

% ----- Binning to obtain low-resolution image --------------------------
imageBlurredCrop = imageHighResBlurred( ...
    1:(lowResSize * zoomFactor), 1:(lowResSize * zoomFactor));

imageLowRes = binning(imageBlurredCrop, ...
    [zoomFactor, zoomFactor], "mean");

fprintf("   - Binning factor: %d x %d\n", zoomFactor, zoomFactor);
fprintf("   - Low-resolution size: %d x %d pixels\n", ...
    size(imageLowRes, 1), size(imageLowRes, 2));

% Optional Gaussian noise in measurement space
imageLowResNoisy = imageLowRes + NOISE_STD * randn(size(imageLowRes));
imageLowResNoisy(imageLowResNoisy < 0) = 0;

if NOISE_STD > 0
    fprintf("   - Added Gaussian noise (std = %.3g)\n", NOISE_STD);
else
    fprintf("   - No measurement noise added (noiseStd = 0)\n");
end

%% 3. Apply ZOOMER for super-resolution reconstruction
fprintf("\nStep 3: Applying ZOOMER super-resolution...\n");

% Reconstruction parameters for zoomer
kPenalty     = 0.005;     % L2 (ridge) penalty on X
lambdaSmooth = 0.0005;    % Smoothness penalty (set > 0 to enforce smoothness)
deltaTol     = 1.0e-10;   % Convergence threshold for CG residual
maxIter      = 1500;      % Maximum number of iterations
sigmaPsfEst  = SIGMA_PSF_SIM;  % PSF sigma used inside zoomer

tic;
[imageReconstructed, numIter, errHist] = zoomer( ...
    imageLowResNoisy, ...
    kPenalty, ...
    lambdaSmooth, ...
    deltaTol, ...
    zoomFactor, ...
    maxIter, ...
    sigmaPsfEst);
elapsedTime = toc;

fprintf("   - Reconstruction completed in %.3f seconds\n", elapsedTime);
fprintf("   - Number of iterations: %d / %d\n", numIter, maxIter);
fprintf("   - Final data-space RMS error: %.6e\n", errHist(end));
fprintf("   - Reconstructed size: %d x %d pixels\n", ...
    size(imageReconstructed, 1), size(imageReconstructed, 2));

%% 4. Display original, low-resolution blurred, and reconstructed images
fprintf("\nStep 4: Displaying results (three images)...\n");

figure("Name", "Zoomer Super-Resolution Demo");
set(gcf, "Position", [100, 100, 1500, 500]);

colormap("hot");

% ----- High-resolution ground truth ------------------------------------
subplot(1, 3, 1);
imagesc(imageHighRes);
axis image off;
title("Ground Truth (High-Res)");

% ----- Low-resolution blurred (measurement) image ----------------------
subplot(1, 3, 2);
imagesc(imageLowResNoisy);
axis image off;
title(sprintf("Low-Res Blurred (zoom %dx)", zoomFactor));

% ----- Reconstructed high-resolution image -----------------------------
subplot(1, 3, 3);
imagesc(imageReconstructed);
axis image off;
title("Superresolved image (clipped neg. entries)");
clim([0 max(imageReconstructed(:))])

fprintf("\n=== Demo completed successfully ===\n");
fprintf("A single figure with three images has been generated.\n\n");
