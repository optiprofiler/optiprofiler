classdef TestProfileBands < matlab.unittest.TestCase
% Rendered profile bands (performance and data profiles with n_runs > 1).
%
% The drawn band of a solver must cover exactly the positive-x-span intervals
% of its stairs-expanded lower/upper arrays and nothing else. Two identical
% runs therefore have a zero-area band and must paint no translucent pixel.
% The legacy single polygon [x; flip(x)], [lower; flip(upper)] has exactly
% that geometric area, but MATLAB tessellates the self-touching polygon into
% triangles and paints broad spurious regions between its extreme corners
% (raster and vector export alike). The benchmark-based cases therefore
% rasterize the saved figures themselves, away from lines, text and legend;
% the helper and single-run cases check geometry and absence directly.
% All fixtures run serially (n_jobs = 1) with a pinned seed.
    properties (Access = private)
        SourceRoot
        OutputRoot
    end
    methods (TestMethodSetup)
        function isolateFixture(testCase)
            folder = fileparts(mfilename('fullpath'));
            [~, root] = fileattrib(fullfile(folder, '..', '..', '..', '..'));
            testCase.SourceRoot = root.Name;
            testCase.OutputRoot = tempname;
            mkdir(testCase.OutputRoot);
            original_path = path;
            original_directory = pwd;
            registry = getenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY');
            rng_state = rng;
            % Teardowns run last-in first-out: the output directory is removed last, after every owned figure has
            % been closed (registered later by the cases), the working directory left it and path, registry and
            % RNG have been restored; this order holds on success, failed assertions and errors alike.
            testCase.addTeardown(@() removeOutput(testCase.OutputRoot));
            testCase.addTeardown(@() rng(rng_state));
            testCase.addTeardown(@() setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', registry));
            testCase.addTeardown(@() path(original_path));
            testCase.addTeardown(@() cd(original_directory));
            addpath(fullfile(testCase.SourceRoot, 'matlab', 'optiprofiler', 'src'));
            % Three deterministic problems of dimension 2, 3 and 4 (tests/fixtures/profilebands): several jumps at
            % pooled duplicate x, the shape on which the legacy band polygon was rendered with spurious area.
            fixture_root = fullfile(folder, '..', 'fixtures', 'profilebands');
            addpath(fixture_root);
            setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', fullfile(testCase.OutputRoot, 'registry.mat'));
            registerProblemLibrary(struct('name', 'band_fixture', 'root', fixture_root, ...
                'select_function', 'band_fixture_select', 'load_function', 'band_fixture_load'));
            cd(testCase.OutputRoot);
        end
    end
    methods (Test)
        function identicalRunsMeanStdPaintNothing(testCase)
            testCase.zeroAreaCase('meanstd');
        end
        function identicalRunsMinMaxPaintNothing(testCase)
            testCase.zeroAreaCase('minmax');
        end
        function singleRunDrawsNoBand(testCase)
            figs = testCase.runBenchmark('single', 1, 'plain', 'meanstd');
            for k = 1:numel(figs)
                [~, patches] = testCase.openBands(figs{k});
                testCase.verifyEmpty(patches, sprintf('%s: a single run must not draw a band', figs{k}));
            end
        end
        function unequalRunsKeepTheirBand(testCase)
            % perturbed_x0 (pinned seed, perturbation level 1 so the start points differ by the order of the
            % problem scale) gives every run its own start point: the deterministic solvers need different amounts
            % of work per run and the across-run bands are broad. Rendered evidence is paired: the same axes are
            % rasterized with the band patches and after deleting them in memory, so painted interior pixels are
            % attributed to the band patches alone.
            [figs, companion] = testCase.runBenchmark('unequal', 2, 'perturbed_x0', 'meanstd');
            figures_with_area = 0;
            for k = 1:numel(figs)
                [ax, patches] = testCase.openBands(figs{k});
                series = companionSeries(companion, figs{k});
                lines = findobj(ax, 'Type', 'line');
                expected_patches = 0; figure_area = 0; box_area = diff(ax.XLim) * diff(ax.YLim);
                for s = 1:numel(series)
                    [expected_vertices, expected_faces] = expectedFaces(series{s});
                    if isempty(expected_faces)
                        continue;
                    end
                    expected_patches = expected_patches + 1;
                    figure_area = figure_area + facesArea(expected_vertices, expected_faces);
                    line = lines(strcmp({lines.DisplayName}, solverName(series{s}.solver_index)));
                    testCase.assertNumElements(line, 1, 'One mean line per solver.');
                    testCase.assertGreaterThanOrEqual(numel(patches), expected_patches, sprintf('%s: band patch of solver %d missing', figs{k}, s));
                    actual = patches(expected_patches);   % patches are created in solver order, only for nonempty bands
                    testCase.verifyEqual(actual.Faces, expected_faces, sprintf('%s solver %d: faces', figs{k}, s));
                    testCase.verifyEqual(actual.Vertices, expected_vertices, sprintf('%s solver %d: vertices', figs{k}, s));
                    testCase.verifyEqual(actual.FaceColor, line.Color, sprintf('%s solver %d: band color', figs{k}, s));
                    testCase.verifyEqual(actual.FaceAlpha, 0.2); testCase.verifyEqual(actual.EdgeColor, 'none');
                    testCase.verifyEqual(actual.HandleVisibility, 'off');
                end
                testCase.verifyNumElements(patches, expected_patches, sprintf('%s: one patch per nonempty band, none otherwise', figs{k}));
                fraction_on = coloredInteriorFraction(ax, fullfile(testCase.OutputRoot, sprintf('unequal-%d-on.png', k)));
                delete(patches);   % in-memory control: the same axes without its band patches
                fraction_off = coloredInteriorFraction(ax, fullfile(testCase.OutputRoot, sprintf('unequal-%d-off.png', k)));
                share = figure_area / box_area;
                testCase.verifyLessThan(fraction_off, 1e-4, sprintf('%s: without band patches nothing may be painted', figs{k}));
                if figure_area > 0
                    figures_with_area = figures_with_area + 1;
                    % The fixture must give broad bands wherever it gives a band (at least 1 % of the axes box), so
                    % that the rendered evidence stays sensitive to an entirely missing band: the patches must paint
                    % at least a quarter of their geometric share (the exported image is the box plus margins,
                    % and pixels near lines, text and legend are excluded) and not grossly beyond it.
                    testCase.verifyGreaterThanOrEqual(share, 0.01, sprintf('%s: fixture band area %.4f of a %.4f box is too thin for rendered evidence', figs{k}, figure_area, box_area));
                    testCase.verifyGreaterThan(fraction_on - fraction_off, 0.25 * share, sprintf('%s: band patches paint %.4f of the image for a geometric share of %.4f', figs{k}, fraction_on - fraction_off, share));
                    testCase.verifyLessThan(fraction_on, 2 * share + 1e-3, sprintf('%s: band painted beyond its geometry', figs{k}));
                else
                    testCase.verifyLessThan(fraction_on, 1e-4, sprintf('%s: no band area, nothing may be painted', figs{k}));
                end
            end
            testCase.verifyGreaterThan(figures_with_area, 0, 'The unequal-run fixture must produce at least one band with positive area.');
        end
        function helperGeometryOnStairsArrays(testCase)
            % Direct geometry of the helper on stairs-expanded arrays (consecutive vertices share x or y).
            % Wholly collapsed band: no face.
            [V, F] = optiprofiler_internal.profileBandFaces([0 1 1 2]', [0 0 0.5 0.5]', [0 0 0.5 0.5]');
            testCase.verifyEmpty(F); testCase.verifyEmpty(V);
            % Mixed band: one positive-area interval among collapsed ones (duplicate x at the jumps).
            [V, F] = optiprofiler_internal.profileBandFaces([0 1 1 2 2 3]', [0 0 0.2 0.2 0.4 0.4]', [0 0 0.8 0.8 0.4 0.4]');
            testCase.verifyEqual(F, [1 2 3 4]); testCase.verifyEqual(V, [1 0.2; 2 0.2; 2 0.8; 1 0.8]);
            % Genuine band built by stairs from sample arrays: zero-height first interval omitted, three faces kept.
            [x, lo] = stairs([0 1 2 3 4]', [0 0.2 0.4 0.6 0.6]'); [~, up] = stairs([0 1 2 3 4]', [0 0.6 0.8 1 1]');
            [V, F] = optiprofiler_internal.profileBandFaces(x, lo, up);
            testCase.verifyEqual(size(F, 1), 3);
            testCase.verifyEqual(V(F(1, :), :), [1 0.2; 2 0.2; 2 0.6; 1 0.6]);
            testCase.verifyEqual(V(F(3, :), :), [3 0.6; 4 0.6; 4 1; 3 1]);
            testCase.verifyEqual(facesArea(V, F), 0.4 + 0.4 + 0.4, 'AbsTol', 1e-12);
            % A strip collapsing at one endpoint keeps only the positive-height intervals.
            [V, F] = optiprofiler_internal.profileBandFaces([0 1 1 2 2 3]', [0 0 0.3 0.3 0.5 0.5]', [0.4 0.4 0.6 0.6 0.5 0.5]');
            testCase.verifyEqual(F, [1 2 3 4; 5 6 7 8]); testCase.verifyEqual(V(5:8, :), [1 0.3; 2 0.3; 2 0.6; 1 0.6]);
            % Nonfinite x: intervals touching it are omitted, nothing is reconnected across the gap.
            [V, F] = optiprofiler_internal.profileBandFaces([0 1 1 NaN 2 2 3]', [0 0 0.2 0.2 0.2 0.4 0.4]', [0.3 0.3 0.5 0.5 0.5 0.9 0.9]');
            testCase.verifyEqual(F, [1 2 3 4; 5 6 7 8]);
            testCase.verifyEqual(V(1:4, :), [0 0; 1 0; 1 0.3; 0 0.3]); testCase.verifyEqual(V(5:8, :), [2 0.4; 3 0.4; 3 0.9; 2 0.9]);
            % Nonfinite y at the RIGHT end vertex of a positive-span interval omits that interval too (both end
            % vertices count): [0,1] touches the NaN at its right end and is dropped, [1,2] is kept.
            [V, F] = optiprofiler_internal.profileBandFaces([0 1 1 2]', [0.2 NaN 0.2 0.2]', [0.6 0.6 0.6 0.6]');
            testCase.verifyEqual(F, [1 2 3 4]); testCase.verifyEqual(V, [1 0.2; 2 0.2; 2 0.6; 1 0.6]);
            % Mismatched lengths are refused, never silently truncated.
            testCase.verifyError(@() optiprofiler_internal.profileBandFaces([0 1 1 2]', [0 0 0.2]', [0 0 0.8 0.8]'), 'OptiProfiler:profileBandFaces:LengthMismatch');
            % Row or column vectors are accepted alike.
            [Vr, Fr] = optiprofiler_internal.profileBandFaces([0 1 1 2], [0 0 0.2 0.2], [0 0 0.8 0.8]);
            testCase.verifyEqual(Fr, [1 2 3 4]); testCase.verifyEqual(Vr, [1 0.2; 2 0.2; 2 0.8; 1 0.8]);
        end
        function adjacentFacesShowNoSeam(testCase)
            % Two faces sharing an edge must paint as one uniform translucent rectangle: the painted run of the middle
            % row must be contiguous across the shared edge (no white crack) and uniform in color (no darker seam).
            f = figure('Visible', 'off'); testCase.addTeardown(@() closeIfValid(f)); ax = axes(f); hold(ax, 'on');
            patch(ax, 'Vertices', [0 0.2; 1 0.2; 1 0.8; 0 0.8; 1 0.2; 2 0.2; 2 0.8; 1 0.8], 'Faces', [1 2 3 4; 5 6 7 8], ...
                'FaceColor', [0 0.447 0.741], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            xlim(ax, [0 2]); ylim(ax, [0 1]); axis(ax, 'off');
            png = fullfile(testCase.OutputRoot, 'seam.png'); exportgraphics(ax, png, 'Resolution', 100);
            image = double(imread(png)); row = squeeze(image(round(size(image, 1) / 2), :, :));
            colored = find(max(row, [], 2) - min(row, [], 2) > 10);
            testCase.assertNotEmpty(colored, 'The band row must be painted.');
            span = colored(1):colored(end);
            testCase.verifyEqual(numel(span), numel(colored), 'The painted run must be contiguous across the shared edge (no white crack).');
            inside = row(span(2:end-1), :);   % the two outermost pixels are anti-aliased against white
            testCase.verifyLessThanOrEqual(max(inside, [], 1) - min(inside, [], 1), [4 4 4], 'The shared edge must not show as a darker seam.');
        end
    end
    methods (Access = private)
        function zeroAreaCase(testCase, errorbar_type)
            [figs, companion] = testCase.runBenchmark(['identical-' errorbar_type], 2, 'plain', errorbar_type);
            for k = 1:numel(figs)
                [ax, patches] = testCase.openBands(figs{k});
                series = companionSeries(companion, figs{k});
                for s = 1:numel(series)
                    [~, expected_faces] = expectedFaces(series{s});
                    testCase.verifyEmpty(expected_faces, sprintf('%s solver %d: identical runs must prepare a zero-area band', figs{k}, s));
                end
                for p = 1:numel(patches)
                    testCase.verifyEqual(facesArea(patches(p).Vertices, patches(p).Faces), 0, 'A zero-area band must not carry faces with area.');
                    testCase.verifyTrue(all(faceSpans(patches(p)) > 0, 'all'), sprintf('%s: every drawn face must have positive x-span and height', figs{k}));
                end
                fraction = coloredInteriorFraction(ax, fullfile(testCase.OutputRoot, sprintf('%s-%d.png', errorbar_type, k)));
                testCase.verifyLessThan(fraction, 1e-4, sprintf('%s: a zero-area band painted %.5f of the rendered axes', figs{k}, fraction));
            end
        end
        function [figs, companion] = runBenchmark(testCase, id, n_runs, feature_name, errorbar_type)
            options = struct('plibs', {{'band_fixture'}}, 'silent', true, 'n_jobs', 1, 'seed', 7, ...
                'feature_name', feature_name, 'n_runs', n_runs, 'score_only', false, 'savepath', testCase.OutputRoot, ...
                'benchmark_id', id, 'draw_hist_plots', 'none', 'max_tol_order', 1, 'max_eval_factor', 60, ...
                'errorbar_type', errorbar_type, 'solver_names', {{'NelderMead', 'CoordWalk'}}, ...
                'report_path', fullfile(testCase.OutputRoot, [id '.json']));
            if strcmp(feature_name, 'perturbed_x0')
                options.perturbation_level = 1;   % start points differ by the order of the problem scale
            end
            benchmark({@nelderMead, @coordWalk}, options);
            report = readJson(fullfile(testCase.OutputRoot, [id '.json']));
            companion = readJson(fullfile(testCase.OutputRoot, report.plot_data.path));
            stamps = dir(fullfile(testCase.OutputRoot, id)); stamps = stamps([stamps.isdir] & ~ismember({stamps.name}, {'.', '..'}));
            testCase.assertNumElements(stamps, 1, 'Exactly one run directory is expected.');
            fig_dir = fullfile(testCase.OutputRoot, id, stamps(1).name, 'test_log', 'profile_figs');
            names = {'perf_hist_1.fig', 'data_hist_1.fig', 'perf_out_1.fig', 'data_out_1.fig'};
            figs = cellfun(@(n) fullfile(fig_dir, n), names, 'UniformOutput', false);
            for k = 1:numel(figs), testCase.assertTrue(isfile(figs{k}), ['Missing figure ' figs{k}]); end
        end
        function [ax, patches] = openBands(testCase, fig_path)
            h = openfig(fig_path, 'new', 'invisible');
            testCase.addTeardown(@() closeIfValid(h));   % closed on success, failure and error alike
            ax = findobj(h, 'Type', 'axes'); ax = ax(1);
            patches = flipud(findall(ax, 'Type', 'patch'));   % creation order: solver order, nonempty bands only
        end
    end
end

function name = solverName(solver_index)
    names = {'NelderMead', 'CoordWalk'}; name = names{solver_index};
end

function series = companionSeries(companion, fig_path)
    [~, name] = fileparts(fig_path); parts = strsplit(name, '_');   % perf|data _ hist|out _ <tolerance index>
    kind = 'performance'; if strcmp(parts{1}, 'data'), kind = 'data'; end
    scope = 'history'; if strcmp(parts{2}, 'out'), scope = 'output'; end
    id = sprintf('profile-%s-%s-%s', parts{3}, scope, kind);
    plots = companion.plots; if isstruct(plots), plots = num2cell(plots); end
    match = plots(cellfun(@(p) strcmp(p.id, id), plots));
    assert(isscalar(match), 'Companion plot %s not found', id);
    series = match{1}.series; if isstruct(series), series = num2cell(series); end
    [~, order] = sort(cellfun(@(s) s.solver_index, series)); series = series(order);
end

function [V, F] = expectedFaces(s)
% Test oracle, written independently of the helper: both end vertices of an interval must be finite.
    x = s.x(:); lo = s.lower(:); up = s.upper(:); V = zeros(0, 2); F = zeros(0, 4);
    for k = 1:numel(x) - 1
        finite = all(isfinite([x(k) x(k+1) lo(k) lo(k+1) up(k) up(k+1)]));
        if finite && x(k+1) > x(k) && up(k) > lo(k)
            base = size(V, 1); V = [V; x(k) lo(k); x(k+1) lo(k); x(k+1) up(k); x(k) up(k)]; F = [F; base + (1:4)]; %#ok<AGROW>
        end
    end
end

function a = facesArea(V, F)
    a = 0;
    for k = 1:size(F, 1)
        q = V(F(k, :), :); a = a + polyarea(q(:, 1), q(:, 2));
    end
end

function spans = faceSpans(p)
    spans = zeros(size(p.Faces, 1), 2);
    for k = 1:size(p.Faces, 1)
        q = p.Vertices(p.Faces(k, :), :); spans(k, :) = [max(q(:, 1)) - min(q(:, 1)), max(q(:, 2)) - min(q(:, 2))];
    end
end

function fraction = coloredInteriorFraction(ax, png)
% Translucent solver color over white (chromatic but light), counted only away from dark pixels: the
% neighbourhoods of mean lines, axes, tick labels and the legend are excluded so that anti-aliased line edges
% cannot masquerade as band interior.
    exportgraphics(ax, png, 'Resolution', 100);
    image = double(imread(png)); mx = max(image, [], 3); mn = min(image, [], 3);
    near_dark = conv2(double(mn < 130), ones(9), 'same') > 0;
    colored = (mx - mn > 10) & (mn > 130) & (mx > 200) & ~near_dark;
    fraction = nnz(colored) / numel(colored);
end

function value = readJson(path)
    fid = fopen(path, 'rb'); bytes = fread(fid, inf, 'uint8=>uint8')'; fclose(fid);
    value = jsondecode(native2unicode(bytes, 'UTF-8'));
end

function x = coordWalk(fun, x0)
    x = x0(:); f = fun(x); step = 0.5;
    for it = 1:20
        improved = false;
        for i = 1:numel(x)
            for sgn = [-1, 1]
                y = x; y(i) = y(i) + sgn * step; fy = fun(y);
                if fy < f, x = y; f = fy; improved = true; end
            end
        end
        if ~improved, step = step / 2; end
    end
end

function x = nelderMead(fun, x0)
    x = fminsearch(fun, x0, optimset('MaxFunEvals', 100, 'Display', 'off'));
end

function closeIfValid(h)
    if isvalid(h), close(h); end
end

function removeOutput(path)
    if isfolder(path), rmdir(path, 's'); end
end
