function [vertices, faces] = profileBandFaces(x, lower, upper)
%PROFILEBANDFACES Rectangular faces of a stairs-expanded across-run band.
%   [VERTICES, FACES] = profileBandFaces(X, LOWER, UPPER) returns one
%   axis-aligned rectangular face per interval of positive x-span and
%   positive height of the band whose stairs-expanded lower and upper
%   boundaries are LOWER and UPPER over X (the outputs of
%   prepareProfilePlotData: consecutive vertices share x or y, and the
%   value carried across an interval [X(k), X(k+1)] with X(k+1) > X(k) is
%   the value at its first vertex k). VERTICES is 4M-by-2, FACES is M-by-4
%   (indices into VERTICES, counterclockwise), both empty when the band has
%   no area. The faces are the rectangles [X(k), X(k+1)] x [LOWER(k),
%   UPPER(k)]. The three arrays must have equal length (stairs outputs of
%   one series). An interval is omitted exactly, without any tolerance, when
%   its x-span or its height is zero (the vertical jumps, including the
%   duplicated x of pooled runs) or when any of the six coordinates of its
%   two end vertices is not finite; nothing is ever connected across an
%   omitted interval.
%
%   Why not one polygon: the legacy fill of [X; flip(X)] with [LOWER;
%   flip(UPPER)] describes the same strip, but when the runs coincide on an
%   interval and differ only at the duplicated x of a jump, the lower and
%   upper chains fold back and forth along the same vertical line. That
%   polygon has zero geometric area there, yet MATLAB tessellates a patch
%   face into triangles and, for such a self-touching polygon, paints broad
%   triangles between its extreme corners, in raster and vector output
%   alike (observed on R2024b and R2026a). Independent rectangles have no
%   self-touching boundary, so the painted set is exactly the strip.
%
%   This function is a presentation helper only: it does not change the
%   prepared arrays, the mean line, the companion plot data or any score.
    x = double(x(:)); lower = double(lower(:)); upper = double(upper(:));
    n = numel(x);
    assert(numel(lower) == n && numel(upper) == n, 'OptiProfiler:profileBandFaces:LengthMismatch', ...
        'X, LOWER and UPPER must be stairs-expanded arrays of equal length (%d, %d and %d given).', n, numel(lower), numel(upper));
    first = 1:n - 1;
    second = 2:n;
    finite = isfinite(x(first)) & isfinite(x(second)) & isfinite(lower(first)) & isfinite(lower(second)) ...
        & isfinite(upper(first)) & isfinite(upper(second));
    keep = find(finite & x(second) > x(first) & upper(first) > lower(first));
    n_faces = numel(keep);
    if n_faces == 0
        vertices = zeros(0, 2);
        faces = zeros(0, 4);
        return;
    end
    % Corner order per face: lower-left, lower-right, upper-right, upper-left.
    corner_x = [x(keep), x(keep + 1), x(keep + 1), x(keep)].';
    corner_y = [lower(keep), lower(keep), upper(keep), upper(keep)].';
    vertices = [corner_x(:), corner_y(:)];
    faces = reshape(1:4 * n_faces, 4, n_faces).';
end
