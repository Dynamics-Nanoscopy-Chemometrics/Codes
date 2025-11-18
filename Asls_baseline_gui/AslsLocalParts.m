% ========================================================================
% AslsLocalParts — Piece-wise ASLS baseline correction (matrix/batch)
% ========================================================================
% This function applies the asymmetric least squares (ASLS) baseline method
% on *selected intervals* for each row of a data matrix and sets the
% baseline to zero elsewhere (with edge blending). Optionally, it performs
% a final ASLS smoothing (p = 0.5) on the obtained baseline.
%
% Method foundation:
%   Eilers, P. H. C., & Boelens, H. F. M. (2005).
%   Baseline correction with asymmetric least squares smoothing.
%   Leiden University Medical Centre Report, 1(1), 5.
%
% ------------------------------
% HOW TO CITE (per Lovelace’s Square)
% ------------------------------
% If you use this code, please cite BOTH:
%
% (1) The original method:
%     Eilers, P. H. C., & Boelens, H. F. M. (2005).
%     Baseline correction with asymmetric least squares smoothing.
%     Leiden University Medical Centre Report, 1(1), 5.
%
% (2) This implementation (example format; update the URL when available):
%     Gómez, Y. R. (2025). AslsLocalParts — piece-wise ASLS (MATLAB).
%     Lovelace’s Square. Version 1.0.0. Last accessed: 2025-11-12.
%     <resource URL on Lovelace’s Square, if available>
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
% ------------------------------
% METADATA
% ------------------------------
% Authors: Yesid Roman Gomez
% Date Created: 2025-10-16
% License: MIT
% Version: 1.0.0
% Implements: Eilers & Boelens (2005) — non-novel implementation
% Dependencies: LocalBaselineOnly, AslsBaseSingle
% Reviewed by Lovelace’s Square team: No
%
% ------------------------------
% WHAT THIS FUNCTION DOES
% ------------------------------
% - Accepts an M×N matrix (rows = spectra/observations, cols = channels).
% - Builds K intervals from startIdx/endIdx (rounded, clamped to [1, N]).
% - For each row, computes a local baseline using LocalBaselineOnly:
%     • ASLS is applied per interval with per-interval ? (lambda) and p.
%     • Outside intervals the baseline is zero (before optional smoothing).
%     • Overlaps are averaged; edges and gaps are blended over 10 samples.
% - Optionally applies a *final* ASLS smoothing to the baseline with p = 0.5
%   and user-specified smoothingLambda. This can make the baseline nonzero
%   outside the intervals due to the smoothing.
% - Returns the per-row baseline and the baseline-corrected data
%   (correctedData = data - baseline), plus a params struct of settings.
%
% ------------------------------
% USAGE
% ------------------------------
% Inputs:
%   data            - [M x N] numeric matrix (rows = spectra, cols = channels)
%   startIdx        - [K x 1] or [1 x K] start indices
%   endIdx          - [K x 1] or [1 x K] end indices
%   lambdas         - [K x 1] or [1 x K] ASLS ? per interval (>= 0)
%   ps              - [K x 1] or [1 x K] ASLS asymmetry p per interval in [0,1]
%   smoothingLambda - scalar; if > 0 applies final ASLS to the baseline with p=0.5
%
% Outputs:
%   correctedData   - [M x N] baseline-corrected data (data - baseline)
%   baseline        - [M x N] estimated baseline, per row
%   params          - struct with effective settings (intervals, ?, p, etc.)
%
% Example:
%   [corr, base, prm] = AslsLocalParts( ...
%       data, [100; 400], [250; 700], [1e5; 1e6], [0.001; 0.01], 0);
%
% ========================================================================

function [correctedData, baseline, params] = AslsLocalParts( ...
        data, startIdx, endIdx, lambdas, ps, smoothingLambda)
% ASLSLOCALPARTS Piece-wise ASLS baseline correction (no GUI)
%
% This function reproduces the behavior of the "Apply Correction" step in
% AslsLocalPartsGui.m, but as a standalone function.
%
% Inputs:
%   data            - [M x N] numeric matrix.
%                     Rows = spectra/observations, Cols = channels.
%   startIdx        - [K x 1] or [1 x K] start indices for intervals.
%   endIdx          - [K x 1] or [1 x K] end indices for intervals.
%   lambdas         - [K x 1] or [1 x K] lambda values per interval.
%   ps              - [K x 1] or [1 x K] p values per interval.
%   smoothingLambda - scalar.
%                     If > 0, apply final ASLS smoothing with p = 0.5.
%                     If <= 0 or empty, no final smoothing (same as
%                     "Enable Smoothing" unchecked).
%
% Outputs:
%   correctedData   - [M x N] baseline-corrected data.
%   baseline        - [M x N] baseline estimated for each row.
%   params          - struct with the effective settings used.
%
% Requirements (same as GUI):
%   - LocalBaselineOnly(y, ints, lams, ps, blendSmooth)
%   - AslsBaseSingle(y, lambda, p)
%
% Notes:
%   - Outside user-defined intervals the baseline is zero; edges are
%     blended over BLEND_SMOOTH samples (as in the GUI via LocalBaselineOnly).
%   - Overlapping intervals are handled/averaged by LocalBaselineOnly.
%
% Authors: Yesid Roman Gomez
% -------------------------------------------------------------------------

    % ---- Constants (matching GUI behavior) ------------------------------
    BLEND_SMOOTH = 10;   % samples for edge blending

    % ---- Input validation -----------------------------------------------
    if nargin < 5
        error('AslsLocalParts:NotEnoughInputs', ...
            'Usage: AslsLocalParts(data, startIdx, endIdx, lambdas, ps, [smoothingLambda])');
    end

    if nargin < 6 || isempty(smoothingLambda)
        smoothingLambda = 0;  % default: no final smoothing
    end

    if ~isnumeric(data) || ndims(data) ~= 2
        error('AslsLocalParts:InvalidData', ...
            'data must be a numeric 2D matrix [M x N].');
    end

    data = double(data);
    [nRows, nCols] = size(data);

    startIdx = startIdx(:);
    endIdx   = endIdx(:);
    lambdas  = lambdas(:);
    ps       = ps(:);

    nInt = numel(startIdx);
    if any([numel(endIdx), numel(lambdas), numel(ps)] ~= nInt)
        error('AslsLocalParts:SizeMismatch', ...
            'startIdx, endIdx, lambdas, ps must have the same length.');
    end

    if nInt == 0
        % No intervals: no baseline, return data unchanged.
        baseline      = zeros(nRows, nCols);
        correctedData = data;
        params = struct( ...
            'localIntervals', [], ...
            'localLambdas', [], ...
            'localPs', [], ...
            'smoothingLambda', 0, ...
            'smoothingEnabled', false, ...
            'blendSmooth', BLEND_SMOOTH);
        return;
    end

    % Build [K x 2] intervals, clamp to [1, nCols], ensure start <= end.
    ints = round([startIdx, endIdx]);
    ints(:, 1) = max(ints(:, 1), 1);
    ints(:, 2) = min(ints(:, 2), nCols);

    bad = ints(:, 1) > ints(:, 2);
    if any(bad)
        error('AslsLocalParts:InvalidIntervals', ...
            'Each interval must satisfy 1 <= start <= end <= N.');
    end

    doSmooth = smoothingLambda > 0;

    % ---- Allocate outputs -----------------------------------------------
    baseline = zeros(nRows, nCols);

    % ---- Row-wise baseline computation (same logic as GUI finalBaseline) -
    runFcn = @(yRow) localFinalBaseline( ...
        yRow, ints, lambdas, ps, BLEND_SMOOTH, doSmooth, smoothingLambda);

    % Try parallel if available; fall back to for-loop.
    try
        parfor r = 1:nRows 
            baseline(r, :) = runFcn(data(r, :));
        end
    catch
        for r = 1:nRows
            baseline(r, :) = runFcn(data(r, :));
        end
    end

    % ---- Corrected data --------------------------------------------------
    correctedData = data - baseline;

    % ---- Exported parameters (mirrors GUI "Export Settings") ------------
    params = struct( ...
        'localIntervals',   ints, ...
        'localLambdas',     lambdas, ...
        'localPs',          ps, ...
        'smoothingEnabled', doSmooth, ...
        'smoothingLambda',  doSmooth * smoothingLambda, ...
        'blendSmooth',      BLEND_SMOOTH);

end

% =====================================================================
% Local helper: replicate GUI finalBaseline() behavior for one row
% =====================================================================
function b = localFinalBaseline(yRow, ints, lams, ps, blendSmooth, ...
        doSmooth, smoothLambda)

    yCol = yRow(:);    % ensure column vector

    % Local piece-wise ASLS baseline (outside intervals zero, blended edges)
    b = LocalBaselineOnly(yCol, ints, lams, ps, blendSmooth);

    % Optional final ASLS smoothing with fixed p=0.5 (as in GUI)
    if doSmooth
        b = AslsBaseSingle(b, smoothLambda, 0.5);
    end

    b = reshape(b, 1, []);   % return as row vector
end
