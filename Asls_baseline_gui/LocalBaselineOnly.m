% ========================================================================
% LocalBaselineOnly — Piece-wise ASLS baseline with zero baseline outside
% ========================================================================
% This code composes a baseline by applying the asymmetric least squares
% (ASLS) method independently on specified intervals and setting the
% baseline to zero elsewhere. It is built on the ASLS approach proposed by:
%
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
%     Gómez, Y. R. (2025). LocalBaselineOnly — piece-wise ASLS (MATLAB).
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
% Dependency: AslsBaseSingle (ASLS solver)
% Reviewed by Lovelace’s Square team: No
%
% ------------------------------
% WHAT THIS FUNCTION DOES
% ------------------------------
% - Runs ASLS independently on each user-specified interval [start, end].
% - Outside those intervals, the baseline is exactly zero.
% - If intervals overlap, their local baselines are averaged where they overlap.
% - To avoid step artifacts at interval boundaries and gaps, linear
%   blending/tapering is applied over 'blendSmooth' samples:
%     • Leading/trailing regions are tapered to/from zero.
%     • Gaps between consecutive intervals are cross-faded.
%
% ------------------------------
% USAGE
% ------------------------------
% Args:
%   y (nx1 double): signal (column or row vector accepted)
%   ints (kx2 double): interval bounds [start, end] (1-based, inclusive)
%   lams (kx1 double): ASLS smoothness ? per interval (>= 0)
%   ps (kx1 double): ASLS asymmetry p per interval, in [0, 1]
%   blendSmooth (scalar): number of samples for edge/gap tapering (>= 0)
%
% Returns:
%   z (nx1 double): piece-wise baseline with blending at edges/gaps
%
% Example:
%   z = LocalBaselineOnly(y, [1,200; 400,700], [1e5;1e6], [0.001;0.01], 10);
% ========================================================================

function z = LocalBaselineOnly(y, ints, lams, ps, blendSmooth)
% LOCALBASELINEONLY Piece-wise ASLS with zero baseline outside intervals

    y = y(:);
    n = numel(y);
    z = zeros(n, 1);

    if isempty(ints)
        return;
    end

    sumB = zeros(n, 1);
    cnt  = zeros(n, 1);

    % Run ASLS in each interval and accumulate
    for k = 1:size(ints, 1)
        s = max(1, ints(k, 1));
        e = min(n, ints(k, 2));
        if s > e
            tmp = s; s = e; e = tmp;
        end

        local = AslsBaseSingle(y(s:e), lams(k), ps(k));
        sumB(s:e) = sumB(s:e) + local;
        cnt(s:e)  = cnt(s:e)  + 1;
    end

    hasVals = cnt > 0;
    z(hasVals) = sumB(hasVals) ./ cnt(hasVals);

    if blendSmooth <= 0
        return;
    end

    % Blend edges and gaps to avoid steps
    [~, ord] = sort(ints(:, 1));
    ints = ints(ord, :);

    blendEdge(1, max(1, ints(1, 1)) - 1, -1);
    for k = 1:size(ints, 1) - 1
        blendGap(ints(k, 2), ints(k + 1, 1));
    end
    blendEdge(min(n, ints(end, 2)) + 1, n, 1);

    % --------------------- nested helpers --------------------------------
    function blendEdge(a, b, dirSign)
        len = min(blendSmooth, abs(b - a + 1));
        if len <= 0, return; end

        if dirSign < 0
            idx = a:(a + len - 1);
            w   = linspace(1, 0, len)';
            z(idx) = w .* z(idx);
        else
            idx = (b - len + 1):b;
            w   = linspace(0, 1, len)';
            z(idx) = w .* z(idx);
        end
    end

    function blendGap(e1, s2)
        if s2 <= e1, return; end
        len = min(blendSmooth, s2 - e1 - 1);
        if len <= 0, return; end

        i1 = (e1 + 1):(e1 + len);
        i2 = (s2 - len):(s2 - 1);

        w = linspace(0, 1, len)';
        z(i1) = w .* z(i1);
        z(i2) = (1 - w) .* z(i2);
    end
end

