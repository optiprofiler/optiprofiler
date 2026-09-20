classdef TestAffineStructure < matlab.unittest.TestCase
% Affine structure safeguard: one decision for bounds and linear constraints.
%
% A feature that changes variables by x = A * y + b (linearly_transformed and
% custom with mod_affine) must hand the solver the SAME problem in the new
% coordinates. There are two representations of the bounds:
%
% - the diagonal shortcut, valid when A is diagonal: the bounds stay bounds,
%   scaled by diag(inv);
% - the generic representation, valid for every invertible A: the solver's
%   bounds become infinite and every finite bound becomes a linear row of A.
%
% The defect these tests pin: the bounds were classified by isdiag(inv) and the
% linear rows by isdiag(A), both exact. A diagonal A whose supplied inverse
% carries any off-diagonal entry, even at roundoff level, made the bounds
% infinite (inv is "not diagonal") while no bound row was added (A "is
% diagonal"). The bounds vanished from the posed problem, silently, while the
% truth channel went on scoring the original bounds.
%
% What is asserted here is the structure the solver receives (xl, xu, aub, bub,
% aeq, beq), not merely that construction returns, and that a point satisfies
% that structure exactly when the truth channel calls it feasible.
%
% The follow-up these tests also pin. The feasible set xl <= A * y + b <= xu is
% a box exactly when A is diagonal, so the decision reads A exactly: an entry
% of 4e-16 next to bounds of 1e16 moves the set by 4, and a tolerance on the
% entries called it roundoff. The inverse cannot change the set and decides
% nothing, except that the shortcut scales by diag(inv) and therefore requires
% it to be the reciprocal of diag(A) to roundoff. One transformation is
% produced per problem and seed, so that user code with a state cannot hand the
% bounds one map and the linear constraints another. A supplied inverse has to
% invert A from both sides, measured so that the units of neither set of
% variables matter. Every finite quantity that is transported (shifted bounds,
% right-hand sides, composed rows, the initial point) raises if it leaves the
% floating-point range, instead of being posed as "no constraint". And the
% data of a problem are double precision numbers whatever class they were
% given in, so that none of this is computed in saturating integer arithmetic.
%
% The second follow-up. The range has two ends: a nonzero bound that falls
% below the smallest normal number, or a transported entry whose terms all do,
% raises as well (1e-200 * [1e-200, 2e-200] was posed as the point [0, 0]). The
% initial point is verified where it is used, because no tolerance on the
% matrices bounds an error at a point: inv * (x0 - b) has to be mapped back to
% x0 by A to the rounding of that evaluation, or the equation is solved from A,
% or construction raises. The derivative methods of a single-feature problem
% follow the change of variables by the chain rule (they returned the
% derivatives of the original callbacks at the solver's point), and a
% composition still provides none. The consistency of a supplied inverse allows
% the rounding of its products and nothing more. And three decisions are pinned
% as they are: what the deprecated conveniences keep (nothing), what mod_bounds
% replaces (the bounds, not the framework's rows), and how an integer beyond
% 2^53 is converted (rounded, as every number is).

    properties (Constant)
        % A diagonal change of variables with a negative entry (bounds must
        % swap) and three scales; every entry and reciprocal is a power of two.
        D = [2; -4; 0.5]
        B = [0.5; -0.5; 1]
        Dense = [2, 0, 1; 1, 1, 0; 0, 0, 1]
        Reference = struct('merit', 0, 'kind', 'lower_bound', 'source', 'sum of squares', ...
            'mapping', 'feasible_objective/1')
        ProblemNames = {'bounded', 'linear', 'fixed'}
    end

    methods (Static)
        function f = sphere(x)
            f = sum(x(:).^2);
        end

        function problem = makeProblem(name)
            s = struct('fun', @TestAffineStructure.sphere, 'x0', [0.5; 0.5; 0.5], ...
                'xl', [-1; -2; -Inf], 'xu', [3; Inf; 4], 'reference', TestAffineStructure.Reference);
            switch name
                case 'bounded'   % finite and infinite bounds, no linear constraints
                case 'linear'    % one linear inequality and one linear equality
                    s.aub = [1, 1, 0]; s.bub = 5; s.aeq = [1, 0, -1]; s.beq = 0.25;
                case 'fixed'     % the first variable is fixed: the generic path must pose an equality row
                    s.x0 = [1; 0.5; 0.5]; s.xl(1) = 1; s.xu(1) = 1; s.aub = [1, 1, 0]; s.bub = 5;
            end
            problem = Problem(s);
        end

        % ------------------------------------------------------------ transforms
        function [A, b, inv] = exactDiagonal(~, ~)
            A = diag(TestAffineStructure.D); b = TestAffineStructure.B; inv = diag(1 ./ TestAffineStructure.D);
        end
        function [A, b, inv] = roundoffInverse(~, ~)
            % Diagonal A; the inverse carries off-diagonal entries of roundoff
            % size (a fraction of one unit in the last place of its diagonal).
            [A, b, inv] = TestAffineStructure.exactDiagonal();
            inv(1, 2) = 1e-17; inv(3, 1) = -0.5 * eps * 0.5;
        end
        function [A, b, inv] = roundoffMatrix(~, ~)
            % The mirror case: the inverse is exactly diagonal and A carries the roundoff.
            [A, b, inv] = TestAffineStructure.exactDiagonal();
            A(2, 1) = 0.5 * eps * 2; A(1, 3) = -1e-17;
        end
        function [A, b, inv] = sloppyInverse(~, ~)
            % Accurate to about 1e-12: far above roundoff and far below the
            % 1e-8 at which the pair would be called inconsistent.
            [A, b, inv] = TestAffineStructure.exactDiagonal();
            inv(1, 2) = 1e-12;
        end
        function [A, b, inv] = sloppyMatrix(~, ~)
            % The mirror case: the inverse is exactly diagonal and A carries a
            % coupling of 1e-12, which is posed.
            [A, b, inv] = TestAffineStructure.exactDiagonal();
            A(2, 1) = 1e-12;
        end
        function [A, b, inv] = couplingHiddenByScale(~, ~)
            % The inverse (diagonal 0.5, -0.25, 2) couples the rows of scales
            % 0.25 and 2 by 2 * eps. That is above n * eps times the smaller of
            % the two scales (0.75 * eps), although below n * eps times the
            % largest entry of the matrix (6 * eps): a tolerance scaled by a
            % norm would overlook it.
            [A, b, inv] = TestAffineStructure.exactDiagonal();
            inv(2, 3) = 2 * eps;
        end
        function [A, b, inv] = dense(~, ~)
            A = TestAffineStructure.Dense; b = TestAffineStructure.B; inv = A \ eye(3);
        end
        function [A, b, inv] = singular(~, ~)
            % A singular A with its pseudo-inverse, the most plausible stand-in
            % for an inverse that does not exist.
            A = [1, 2, 0; 2, 4, 0; 0, 0, 1]; b = zeros(3, 1); inv = [0.04, 0.08, 0; 0.08, 0.16, 0; 0, 0, 1];
        end
        function [A, b, inv] = identityAsInverse(~, ~)
            A = TestAffineStructure.Dense; b = TestAffineStructure.B; inv = eye(3);
        end
        function [A, b, inv] = materiallyWrongInverse(~, ~)
            [A, b, inv] = TestAffineStructure.exactDiagonal(); inv(1, 2) = 1e-3;
        end
        function [A, b, inv] = nanInMatrix(~, ~)
            [A, b, inv] = TestAffineStructure.exactDiagonal(); A(1, 2) = NaN;
        end
        function [A, b, inv] = infInInverse(~, ~)
            [A, b, inv] = TestAffineStructure.exactDiagonal(); inv(2, 3) = Inf;
        end
        function [A, b, inv] = nanInShift(~, ~)
            [A, b, inv] = TestAffineStructure.exactDiagonal(); b(2) = NaN;
        end
        function [A, b, inv] = wrongMatrixShape(~, ~)
            A = eye(2); b = zeros(3, 1); inv = eye(3);
        end
        function [A, b, inv] = wrongInverseShape(~, ~)
            A = eye(3); b = zeros(3, 1); inv = eye(2);
        end
        function [A, b, inv] = wrongShiftSize(~, ~)
            A = eye(3); b = zeros(2, 1); inv = eye(3);
        end
        function [A, b, inv] = complexMatrix(~, ~)
            [A, b, inv] = TestAffineStructure.exactDiagonal(); A = complex(A, 0); A(1, 1) = 2 + 1e-3i;
        end
        function [A, b, inv] = scaledRotation(ratio)
            % A rotation after scaling one NEW variable by RATIO, with its
            % inverse (consistent to roundoff).
            rotation = [2, -2, 1; 1, 2, 2; 2, 1, -2] / 3;
            scale = [1, ratio, 1];
            A = rotation .* scale; b = zeros(3, 1); inv = (rotation ./ scale)';
        end
        function [A, b, inv] = numericallySingular(~, ~)
            % Consistent to roundoff, so the residual test cannot see it, and
            % singular to working precision: no scaling of the original
            % variables repairs 1e17.
            [A, b, inv] = TestAffineStructure.scaledRotation(1e17);
        end
        function [A, b, inv] = badlyScaledRotation(~, ~)
            % The same with 1e13: badly scaled, invertible, and to be accepted.
            [A, b, inv] = TestAffineStructure.scaledRotation(1e13);
        end
        function [A, b, inv] = extremeDiagonal(~, ~)
            % An exact diagonal scaling over 300 orders of magnitude: its usual
            % condition number is 1e300, its posed problem exact.
            d = [1e-150; -4; 1e150];
            A = diag(d); b = TestAffineStructure.B; inv = diag(1 ./ d);
        end
        function [A, b, inv] = overflowingScale(~, ~)
            % Exactly consistent and perfectly conditioned, yet the scaled
            % bound 2^1023 * 3 overflows: a finite bound would become infinite.
            A = 2^-1023 * eye(3); b = zeros(3, 1); inv = 2^1023 * eye(3);
        end
        function [A, b, inv] = roundingLevelCoupling(~, ~)
            % A differs from the identity by one entry of 4e-16, which is below
            % n * eps, and the inverse is exact. With bounds of 1e16 that entry
            % moves the feasible set by 4: a coupling, however small it looks.
            A = [1, 4e-16; 0, 1]; b = zeros(2, 1); inv = [1, -4e-16; 0, 1];
        end
        function [A, b, inv] = inaccurateDiagonalInverse(~, ~)
            % Exactly diagonal, and consistent to 1e-9, which the contract
            % accepts (1e-8). Bounds scaled by this diag(inv) would be off by
            % 1e-9 of their size.
            [A, b] = TestAffineStructure.exactDiagonal(); inv = diag((1 + 1e-9) ./ TestAffineStructure.D);
        end
        function [A, b, inv] = oneSidedInverse(~, ~)
            % A * inv is the identity to 1e-16, and inv * A misses it by 1:
            % inv * (A * y) is not y.
            A = diag([1e-8, 1e8, 1]); b = zeros(3, 1); inv = diag([1e8, 1e-8, 1]); inv(1, 2) = 1e-8;
        end
        function [A, b, inv] = oneSidedInverseMirror(~, ~)
            % The same with the large and the small scale exchanged, so that
            % the stray entry sits below the diagonal.
            A = diag([1e8, 1e-8, 1]); b = zeros(3, 1); inv = diag([1e-8, 1e8, 1]); inv(2, 1) = 1e-8;
        end
        function rotation = genericRotation()
            % A rotation that is orthogonal to roundoff only. The entries of
            % the rotation of scaledRotation are t and 2 * t, so that its
            % products cancel exactly in plain floating-point arithmetic and
            % not in fused arithmetic: whether the base accepted it when scaled
            % depended on the kernel that multiplied it.
            c1 = cos(0.3); s1 = sin(0.3); c2 = cos(0.5); s2 = sin(0.5);
            rotation = [c1, -s1, 0; s1, c1, 0; 0, 0, 1] * [1, 0, 0; 0, c2, -s2; 0, s2, c2];
        end
        function [A, b, inv] = rowScaledRotation(~, ~)
            % The mirror of badlyScaledRotation: one ORIGINAL variable is in
            % units of 1e13 (a row of A). Consistent to roundoff and well
            % conditioned in the sense that is checked: to be accepted as well.
            rotation = TestAffineStructure.genericRotation();
            scale = [1; 1e13; 1];
            A = rotation .* scale; b = zeros(3, 1); inv = rotation' ./ scale';
        end
        function [A, b, inv] = columnScaledRotation(~, ~)
            % The same rotation with one NEW variable in units of 1e13 (a column of A).
            rotation = TestAffineStructure.genericRotation();
            scale = [1, 1e13, 1];
            A = rotation .* scale; b = zeros(3, 1); inv = (rotation ./ scale)';
        end
        function [A, b, inv] = shiftedOutOfRange(~, ~)
            A = eye(3); b = -1e308 * ones(3, 1); inv = eye(3);
        end
        function [A, b, inv] = denseShiftedOutOfRange(~, ~)
            A = TestAffineStructure.Dense; b = -1e308 * ones(3, 1); inv = A \ eye(3);
        end
        function [A, b, inv] = shiftedTheOtherWay(~, ~)
            A = eye(3); b = 1e308 * ones(3, 1); inv = eye(3);
        end
        function [A, b, inv] = hugeScaling(~, ~)
            % Exactly consistent and perfectly conditioned; a coefficient of
            % 1e200 times 1e200 overflows.
            A = diag([1e200, 1, 1]); b = zeros(3, 1); inv = diag([1e-200, 1, 1]);
        end

        function transform = powerScaling(exponent)
            % x = 2^exponent * y in the first variable, with its exact inverse
            % (powers of two: nothing is rounded).
            transform = @(s, p) deal(diag([2^exponent, 1, 1]), zeros(3, 1), diag([2^-exponent, 1, 1]));
        end
        function [A, b, inv] = denseWithShift(~, ~)
            A = TestAffineStructure.Dense; b = TestAffineStructure.B; inv = A \ eye(3);
        end
        function problem = smoothProblem(derivatives)
            % A smooth problem with every derivative (or with none of them).
            s = struct('fun', @(x) sum(x.^4) + x(1) * x(2) + sin(x(3)), 'x0', [0.3; -0.4; 0.5], ...
                'xl', [-2; -2; -2], 'xu', [2; 2; 2], ...
                'cub', @(x) [x(1)^2 + x(2) * x(3) - 1; exp(x(1)) - x(3)], 'ceq', @(x) x(1) * x(2) * x(3) - 0.5);
            if derivatives
                s.grad = @(x) 4 * x(:).^3 + [x(2); x(1); cos(x(3))];
                s.hess = @(x) diag(12 * x(:).^2) + [0, 1, 0; 1, 0, 0; 0, 0, -sin(x(3))];
                s.jcub = @(x) [2 * x(1), x(3), x(2); exp(x(1)), 0, -1];
                s.hcub = @(x) {[2, 0, 0; 0, 0, 1; 0, 1, 0]; [exp(x(1)), 0, 0; 0, 0, 0; 0, 0, 0]};
                s.jceq = @(x) [x(2) * x(3), x(1) * x(3), x(1) * x(2)];
                s.hceq = @(x) {[0, x(3), x(2); x(3), 0, x(1); x(2), x(1), 0]};
            end
            problem = Problem(s);
        end
        function features = variableChanges()
            % Every single feature that changes the variables.
            features = {Feature('custom', struct('mod_affine', @TestAffineStructure.denseWithShift)), ...
                Feature('custom', struct('mod_affine', @TestAffineStructure.exactDiagonal)), ...
                Feature('linearly_transformed', struct('rotated', true, 'condition_factor', 4)), ...
                Feature('permuted')};
        end
        function [A, b, inv] = transformationOf(feature, seed, problem)
            % The transformation that a featured problem of this single feature
            % and seed is built with (a kernel of the same stage and seed).
            stage = feature.stages{1};
            kernel = optiprofiler_internal.FeatureKernel(stage.name, stage.options);
            [A, b, inv] = kernel.modifier_affine(seed, problem);
        end
        function row = rowOf(matrix, i)
            row = matrix(i, :);
        end
        function J = centralDifferences(f, y)
            % One row per output, one column per variable.
            h = 1e-5; J = [];
            for k = 1:numel(y)
                e = zeros(numel(y), 1); e(k) = h;
                column = (f(y + e) - f(y - e)) / (2 * h);
                J(:, k) = column(:); %#ok<AGROW>
            end
        end

        function stage = customStage(transform)
            stage = struct('name', 'custom', 'options', struct('mod_affine', transform));
        end

        function stages = pipeline(transform, composed)
            % One custom stage with mod_affine, followed by 'noisy' in a composition.
            stages = {TestAffineStructure.customStage(transform)};
            if composed
                stages = [stages, {'noisy'}];
            end
        end

        % ------------------------------------------------ what the solver receives
        function posed = posedOf(featured)
            posed = struct('xl', featured.xl, 'xu', featured.xu, 'aub', featured.aub, 'bub', featured.bub, ...
                'aeq', featured.aeq, 'beq', featured.beq);
        end

        function posed = posedByKernel(kernel, seed, problem)
            [xl, xu] = kernel.modifier_bounds(seed, problem);
            [aub, bub] = kernel.modifier_linear_ub(seed, problem);
            [aeq, beq] = kernel.modifier_linear_eq(seed, problem);
            posed = struct('xl', xl, 'xu', xu, 'aub', aub, 'bub', bub, 'aeq', aeq, 'beq', beq);
        end

        function expected = expectedDiagonal(problem, A, b, inv)
            % The diagonal shortcut: scaled bounds (swapped where the scale is
            % negative) and linear maps composed with A.
            scale = diag(inv);
            lower = scale .* (problem.xl - b); upper = scale .* (problem.xu - b);
            expected = struct('xl', min(lower, upper), 'xu', max(lower, upper), ...
                'aub', problem.aub * A, 'bub', problem.bub - problem.aub * b, ...
                'aeq', problem.aeq * A, 'beq', problem.beq - problem.aeq * b);
        end

        function expected = expectedGeneric(problem, A, b)
            % The generic representation: no solver bounds; every finite bound
            % is a row of A (an equality row if the variable is fixed).
            fixed = problem.xl == problem.xu;
            upper = isfinite(problem.xu) & ~fixed; lower = isfinite(problem.xl) & ~fixed;
            expected = struct('xl', -Inf(problem.n, 1), 'xu', Inf(problem.n, 1), ...
                'aub', [A(upper, :); -A(lower, :); problem.aub * A], ...
                'bub', [problem.xu(upper) - b(upper); -(problem.xl(lower) - b(lower)); problem.bub - problem.aub * b], ...
                'aeq', [A(fixed, :); problem.aeq * A], ...
                'beq', [problem.xu(fixed) - b(fixed); problem.beq - problem.aeq * b]);
        end

        function [M, c] = mapOf(featured)
            % The affine map from solver to original coordinates, recovered
            % from n + 1 evaluations of the hidden coordinate map.
            n = featured.n;
            c = featured.toOriginalCoordinates(zeros(n, 1));
            M = zeros(n); basis = eye(n);
            for j = 1:n
                M(:, j) = featured.toOriginalCoordinates(basis(:, j)) - c;
            end
        end

        function violation = posedViolation(posed, y)
            % Largest violation of the structure handed to the solver.
            violation = 0;
            lower = isfinite(posed.xl); upper = isfinite(posed.xu);
            if any(lower), violation = max(violation, max(posed.xl(lower) - y(lower))); end
            if any(upper), violation = max(violation, max(y(upper) - posed.xu(upper))); end
            if ~isempty(posed.aub), violation = max(violation, max(posed.aub * y - posed.bub)); end
            if ~isempty(posed.aeq), violation = max(violation, max(abs(posed.aeq * y - posed.beq))); end
        end
    end

    methods
        function verifyStructure(testCase, posed, expected, label)
            names = fieldnames(expected);
            for k = 1:numel(names)
                testCase.verifyEqual(posed.(names{k}), expected.(names{k}), [label, ' ', names{k}]);
            end
        end

        function verifyNoBoundIsLostOrDoubled(testCase, posed, problem, label)
            % Every finite original bound is posed exactly once: as a solver
            % bound or as a linear row, never neither and never both.
            finite = nnz(isfinite(problem.xl)) + nnz(isfinite(problem.xu));
            as_bounds = nnz(isfinite(posed.xl)) + nnz(isfinite(posed.xu));
            as_rows = (size(posed.aub, 1) - problem.m_linear_ub) + 2 * (size(posed.aeq, 1) - problem.m_linear_eq);
            testCase.verifyEqual(as_bounds + as_rows, finite, [label, ' every finite bound is posed exactly once']);
            testCase.verifyTrue(as_bounds == 0 || as_bounds == finite, [label, ' all solver bounds or all linear rows']);
        end

        function verifyPosedIsScored(testCase, posed, A, b, truth, problem, label)
            % A point satisfies the solver-facing structure exactly when the
            % truth (the original problem at A * y + b) calls it feasible. A
            % dropped bound breaks this at once: far outside the original box
            % the posed structure reports no violation and the truth a large one.
            stream = RandStream('mt19937ar', 'Seed', 0);
            points = -6 + 14 * rand(stream, problem.n, 400);
            feasible = 0; infeasible = 0; disagreements = 0;
            fixed = problem.xl == problem.xu;
            for j = 1:size(points, 2)
                x = points(:, j);
                if problem.m_linear_eq > 0   % sample on the equality so that feasible points exist
                    x = x - problem.aeq' * ((problem.aeq * problem.aeq') \ (problem.aeq * x - problem.beq));
                end
                x(fixed) = problem.xl(fixed);
                y = A \ (x - b);
                scored = truth(y); seen = TestAffineStructure.posedViolation(posed, y);
                if (scored > 1e-9 && scored < 1e-3) || (seen > 1e-9 && seen < 1e-3)
                    continue   % too close to a boundary to classify in floating point
                end
                disagreements = disagreements + ((scored <= 1e-9) ~= (seen <= 1e-9));
                feasible = feasible + (scored <= 1e-9); infeasible = infeasible + (scored > 1e-9);
            end
            testCase.verifyEqual(disagreements, 0, [label, ' posed feasibility differs from scored feasibility']);
            testCase.verifyGreaterThan(feasible, 20, label);   % the sample covers both sides
            testCase.verifyGreaterThan(infeasible, 20, label);
        end

        function verifyFeatured(testCase, featured, problem, expected, label)
            posed = TestAffineStructure.posedOf(featured);
            if ~isempty(expected), testCase.verifyStructure(posed, expected, label); end
            testCase.verifyNoBoundIsLostOrDoubled(posed, problem, label);
            [M, c] = TestAffineStructure.mapOf(featured);
            testCase.verifyPosedIsScored(posed, M, c, @(y) scoredViolation(featured, y), problem, label);
        end

        function useAffineFixtures(testCase)
            % PerturbedAffineKernel and StatefulAffine, for this test only.
            old_path = path;
            testCase.addTeardown(@() path(old_path));
            addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'fixtures', 'affine'));
        end
    end

    methods (Test)
        % ------------------------------------------------------ custom with mod_affine

        function exactDiagonalKeepsBoundsAsBounds(testCase)
            [A, b, inv] = TestAffineStructure.exactDiagonal();
            for k = 1:numel(TestAffineStructure.ProblemNames)
                name = TestAffineStructure.ProblemNames{k};
                problem = TestAffineStructure.makeProblem(name);
                featured = FeaturedProblem(problem, Feature('custom', struct('mod_affine', @TestAffineStructure.exactDiagonal)), 10, 3);
                testCase.verifyFeatured(featured, problem, TestAffineStructure.expectedDiagonal(problem, A, b, inv), name);
                % The negative scale swapped the second pair of bounds: (-2, Inf) became (-Inf, 0.375).
                testCase.verifyEqual(featured.xu(2:3), [0.375; 6], name);
                testCase.verifyEqual(featured.xl(2:3), [-Inf; -Inf], name);
                testCase.verifyEqual(A * featured.x0 + b, problem.x0, 'AbsTol', 1e-15, name);
            end
            testCase.verifyEqual(featured.xl(1), 0.25);  % the fixed variable: 0.5 * (1 - 0.5)
        end

        function strayEntriesOfTheInverseKeepBoundsAsBounds(testCase)
            % The original defect: with a stray entry in the inverse the bounds
            % became infinite and no bound row was added, so the bounds were
            % gone. A is exactly diagonal in all three, so the feasible set is
            % a box whatever the inverse carries off its diagonal, at roundoff
            % level or above it: the inverse cannot change the set and decides
            % nothing.
            transforms = {@TestAffineStructure.roundoffInverse, @TestAffineStructure.sloppyInverse, ...
                @TestAffineStructure.couplingHiddenByScale};
            for t = 1:numel(transforms)
                [A, b, inv] = transforms{t}();
                for k = 1:numel(TestAffineStructure.ProblemNames)
                    name = TestAffineStructure.ProblemNames{k};
                    label = sprintf('%s %s', func2str(transforms{t}), name);
                    problem = TestAffineStructure.makeProblem(name);
                    featured = FeaturedProblem(problem, Feature('custom', struct('mod_affine', transforms{t})), 10, 3);
                    testCase.verifyFeatured(featured, problem, TestAffineStructure.expectedDiagonal(problem, A, b, inv), label);
                    % Every finite bound is still a finite solver bound.
                    testCase.verifyEqual(nnz(isfinite(featured.xl)) + nnz(isfinite(featured.xu)), ...
                        nnz(isfinite(problem.xl)) + nnz(isfinite(problem.xu)), label);
                    % The reference is retained, and rightly so: the posed problem is the scored one.
                    testCase.verifyEqual(featured.reference, TestAffineStructure.Reference, label);
                end
            end
        end

        function anyCouplingInTheMatrixTakesTheGenericPath(testCase)
            % An off-diagonal entry of A couples two variables, so the feasible
            % set is no box, and how far it is from one depends on the size of
            % the other variable, which no tolerance on the entry knows. It is
            % posed, exactly, by the representation that needs only A. On the
            % base a coupling in the matrix posed every bound twice.
            transforms = {@TestAffineStructure.roundoffMatrix, @TestAffineStructure.sloppyMatrix};
            for t = 1:numel(transforms)
                [A, b] = transforms{t}();
                for k = 1:numel(TestAffineStructure.ProblemNames)
                    name = TestAffineStructure.ProblemNames{k};
                    label = [func2str(transforms{t}), ' ', name];
                    problem = TestAffineStructure.makeProblem(name);
                    featured = FeaturedProblem(problem, Feature('custom', struct('mod_affine', transforms{t})), 10, 3);
                    testCase.verifyFeatured(featured, problem, TestAffineStructure.expectedGeneric(problem, A, b), label);
                    testCase.verifyEqual(featured.reference, TestAffineStructure.Reference, label);
                end
            end
        end

        function badlyScaledTransformationsKeepWorking(testCase)
            % Scaling is not singularity. The condition number that is checked
            % does not change with the units of the original variables, so an
            % exact diagonal scaling passes whatever its entries, and a
            % rotation of badly scaled variables passes while it is invertible
            % in double precision. Structures are compared exactly; sampling
            % points at these scales would measure roundoff, not structure.
            for k = 1:numel(TestAffineStructure.ProblemNames)
                name = TestAffineStructure.ProblemNames{k};
                problem = TestAffineStructure.makeProblem(name);
                [A, b, inv] = TestAffineStructure.extremeDiagonal();
                featured = FeaturedProblem(problem, Feature('custom', struct('mod_affine', @TestAffineStructure.extremeDiagonal)), 10, 3);
                posed = TestAffineStructure.posedOf(featured);
                testCase.verifyStructure(posed, TestAffineStructure.expectedDiagonal(problem, A, b, inv), ['diagonal ', name]);
                testCase.verifyNoBoundIsLostOrDoubled(posed, problem, ['diagonal ', name]);
                [A, b] = TestAffineStructure.badlyScaledRotation();
                featured = FeaturedProblem(problem, Feature('custom', struct('mod_affine', @TestAffineStructure.badlyScaledRotation)), 10, 3);
                posed = TestAffineStructure.posedOf(featured);
                testCase.verifyStructure(posed, TestAffineStructure.expectedGeneric(problem, A, b), ['rotation ', name]);
                testCase.verifyNoBoundIsLostOrDoubled(posed, problem, ['rotation ', name]);
                % The units of the ORIGINAL variables matter as little as those
                % of the new ones. A * inv misses the identity by 5e-4 here,
                % roundoff times the ratio of the units, and the base called
                % that inconsistent although the condition number it checks is 2.8.
                [A, b, inv] = TestAffineStructure.rowScaledRotation();
                testCase.verifyGreaterThan(norm(A * inv - eye(3), 'fro'), 1e-8 * 3, name);
                featured = FeaturedProblem(problem, Feature('custom', struct('mod_affine', @TestAffineStructure.rowScaledRotation)), 10, 3);
                posed = TestAffineStructure.posedOf(featured);
                testCase.verifyStructure(posed, TestAffineStructure.expectedGeneric(problem, A, b), ['row scaled ', name]);
                testCase.verifyNoBoundIsLostOrDoubled(posed, problem, ['row scaled ', name]);
            end
        end

        function aProblemWithoutVariablesIsStillBuilt(testCase)
            % Empty arrays come in more than one shape (0-by-0, 0-by-1); the
            % validation must not be what breaks the empty problem.
            problem = Problem(struct('fun', @TestAffineStructure.sphere, 'x0', zeros(0, 1)));
            stage = TestAffineStructure.customStage(@(s, p) deal(zeros(0, 0), zeros(0, 1), zeros(0, 0)));
            featured = FeaturedProblem(problem, Feature({stage}), 5, 3);
            testCase.verifyEqual(featured.n, 0);
            testCase.verifyEmpty(featured.xl);
            testCase.verifyEmpty(featured.aub);
            featured = FeaturedProblem(problem, Feature({stage, 'noisy'}), 5, 3);
            testCase.verifyEqual(featured.n, 0);
        end

        function denseTransformKeepsEveryBoundAsALinearRow(testCase)
            [A, b] = TestAffineStructure.dense();
            for k = 1:numel(TestAffineStructure.ProblemNames)
                name = TestAffineStructure.ProblemNames{k};
                problem = TestAffineStructure.makeProblem(name);
                featured = FeaturedProblem(problem, Feature('custom', struct('mod_affine', @TestAffineStructure.dense)), 10, 3);
                testCase.verifyFeatured(featured, problem, TestAffineStructure.expectedGeneric(problem, A, b), name);
                testCase.verifyFalse(any(isfinite(featured.xl)) || any(isfinite(featured.xu)), name);
            end
        end

        function invalidTransformsFailClosed(testCase)
            % No featured problem exists afterwards, so nothing can be posed,
            % scored or claimed. NaN used to pass: norm(...) > tolerance is
            % false for NaN. In a composition the shape and type of a callback
            % output are rejected first, by the custom stage itself.
            invalid = 'MATLAB:Feature:AffineTransformationInvalid';            % A or b is not usable data
            not_invertible = 'MATLAB:Feature:AffineTransformationNotInvertible';  % inv is not the inverse of A
            output = 'MATLAB:Feature:InvalidCustomOutput';
            % Transform, identifier as a single stage, identifier in a
            % composition, and the reason that the single stage must give.
            cases = { ...
                @TestAffineStructure.singular, not_invertible, not_invertible, 'not an identity matrix'; ...
                @TestAffineStructure.identityAsInverse, not_invertible, not_invertible, 'not an identity matrix'; ...
                @TestAffineStructure.materiallyWrongInverse, not_invertible, not_invertible, 'not an identity matrix'; ...
                @TestAffineStructure.nanInMatrix, invalid, invalid, 'must be finite'; ...
                @TestAffineStructure.infInInverse, not_invertible, not_invertible, 'must be finite'; ...
                @TestAffineStructure.nanInShift, invalid, invalid, 'must be finite'; ...
                @TestAffineStructure.wrongMatrixShape, invalid, output, 'real matrix of size 3-by-3'; ...
                @TestAffineStructure.wrongInverseShape, not_invertible, output, 'real matrix of size 3-by-3'; ...
                @TestAffineStructure.wrongShiftSize, invalid, output, 'real vector of size 3'; ...
                @TestAffineStructure.complexMatrix, invalid, output, 'real matrix of size 3-by-3'; ...
                @TestAffineStructure.numericallySingular, not_invertible, not_invertible, 'numerically singular'; ...
                @TestAffineStructure.oneSidedInverse, not_invertible, not_invertible, 'not an identity matrix'; ...
                @TestAffineStructure.oneSidedInverseMirror, not_invertible, not_invertible, 'not an identity matrix'};
            for c = 1:size(cases, 1)
                stage = TestAffineStructure.customStage(cases{c, 1});
                for k = 1:numel(TestAffineStructure.ProblemNames)
                    problem = TestAffineStructure.makeProblem(TestAffineStructure.ProblemNames{k});
                    label = sprintf('%s %s', func2str(cases{c, 1}), TestAffineStructure.ProblemNames{k});
                    testCase.verifyError(@() FeaturedProblem(problem, Feature({stage}), 10, 3), cases{c, 2}, label);
                    testCase.verifyError(@() FeaturedProblem(problem, Feature({stage, 'noisy'}), 10, 3), cases{c, 3}, label);
                    try
                        FeaturedProblem(problem, Feature({stage}), 10, 3);
                    catch exception
                        testCase.verifySubstring(exception.message, cases{c, 4}, label);
                    end
                end
            end
        end

        function finiteBoundNeverBecomesInfiniteByOverflow(testCase)
            % Every check on the matrices passes (the pair is exactly consistent
            % and perfectly conditioned); only the scaled bound overflows.
            % Letting it through would drop the bound without a word.
            problem = TestAffineStructure.makeProblem('bounded');
            feature = Feature('custom', struct('mod_affine', @TestAffineStructure.overflowingScale));
            testCase.verifyError(@() FeaturedProblem(problem, feature, 10, 3), 'MATLAB:Feature:AffineBoundsNotRepresentable');
        end

        function suppliedLinearRowsCannotSilentlyReplaceTheBoundRows(testCase)
            % A supplied linear modifier replaces the linear constraints
            % verbatim (an established rule). Under a dense map the framework
            % poses the bounds as exactly those constraints, so the combination
            % had no place left for the bounds and dropped them. It now raises,
            % and only when a bound is really at stake.
            rows = @(s, p) deal([1, 1, 0], 5);
            equalities = @(s, p) deal([1, 0, -1], 0.25);
            box = @(s, p) deal(-9 * ones(3, 1), 9 * ones(3, 1));
            dense = @TestAffineStructure.dense;
            bounded = TestAffineStructure.makeProblem('bounded');
            fixed = TestAffineStructure.makeProblem('fixed');
            unbounded = Problem(struct('fun', @TestAffineStructure.sphere, 'x0', [0.5; 0.5; 0.5], 'aub', [1, 1, 0], 'bub', 5));
            lost = 'MATLAB:Feature:AffineBoundsNotRepresentable';
            for composed = [false, true]
                label = sprintf('composed=%d', composed);
                tail = {};
                if composed, tail = {'noisy'}; end
                build = @(problem, options) FeaturedProblem(problem, ...
                    Feature([{struct('name', 'custom', 'options', options)}, tail]), 10, 3);
                testCase.verifyError(@() build(bounded, struct('mod_affine', dense, 'mod_linear_ub', rows)), lost, label);
                testCase.verifyError(@() build(fixed, struct('mod_affine', dense, 'mod_linear_eq', equalities)), lost, label);
                % The same decision is read here: a diagonal map whose inverse
                % is not good enough for the shortcut poses the bounds as rows
                % as well.
                testCase.verifyError(@() build(bounded, struct('mod_affine', @TestAffineStructure.inaccurateDiagonalInverse, ...
                    'mod_linear_ub', rows)), lost, label);

                % Nothing fixed, so replacing the equalities loses nothing: the bounds are inequality rows.
                featured = build(bounded, struct('mod_affine', dense, 'mod_linear_eq', equalities));
                expected = TestAffineStructure.expectedGeneric(bounded, TestAffineStructure.Dense, TestAffineStructure.B);
                testCase.verifyEqual(featured.aub, expected.aub, label);
                testCase.verifyEqual(featured.aeq, [1, 0, -1], label);
                % No finite bound, so nothing is at stake.
                featured = build(unbounded, struct('mod_affine', dense, 'mod_linear_ub', rows));
                testCase.verifyEqual(featured.aub, [1, 1, 0], label);
                % The user takes over the bounds as well: both are theirs, verbatim.
                featured = build(bounded, struct('mod_affine', dense, 'mod_linear_ub', rows, 'mod_bounds', box));
                testCase.verifyEqual(featured.xl, -9 * ones(3, 1), label);
                testCase.verifyEqual(featured.aub, [1, 1, 0], label);
                % A diagonal map keeps the bounds as bounds, so supplied rows replace nothing of theirs.
                featured = build(bounded, struct('mod_affine', @TestAffineStructure.exactDiagonal, 'mod_linear_ub', rows));
                testCase.verifyEqual(nnz(isfinite(featured.xl)) + nnz(isfinite(featured.xu)), 4, label);
                testCase.verifyEqual(featured.aub, [1, 1, 0], label);
            end
        end

        function structureDoesNotDependOnWhichModifierIsAskedFirst(testCase)
            % One decision, read by all three modifiers: asking in any order,
            % any number of times, gives the same representation.
            problem = TestAffineStructure.makeProblem('linear');
            transforms = {@TestAffineStructure.exactDiagonal, @TestAffineStructure.roundoffInverse, ...
                @TestAffineStructure.roundoffMatrix, @TestAffineStructure.sloppyInverse, @TestAffineStructure.dense};
            for t = 1:numel(transforms)
                kernel = optiprofiler_internal.FeatureKernel('custom', struct('mod_affine', transforms{t}));
                [aeq1, beq1] = kernel.modifier_linear_eq(3, problem);
                [aub1, bub1] = kernel.modifier_linear_ub(3, problem);
                [xl1, xu1] = kernel.modifier_bounds(3, problem);
                posed = TestAffineStructure.posedByKernel(kernel, 3, problem);
                testCase.verifyEqual({xl1, xu1, aub1, bub1, aeq1, beq1}, ...
                    {posed.xl, posed.xu, posed.aub, posed.bub, posed.aeq, posed.beq}, func2str(transforms{t}));
                bounds_kept = any(isfinite(posed.xl)) || any(isfinite(posed.xu));
                rows_added = size(posed.aub, 1) > problem.m_linear_ub;
                testCase.verifyNotEqual(bounds_kept, rows_added, func2str(transforms{t}));
            end
        end

        % -------------------------------------------------------- linearly_transformed

        function unrotatedIsTheExactDiagonalCaseAndRotatedTheDenseCase(testCase)
            for k = 1:numel(TestAffineStructure.ProblemNames)
                name = TestAffineStructure.ProblemNames{k};
                problem = TestAffineStructure.makeProblem(name);
                options = struct('rotated', false, 'condition_factor', 4);
                featured = FeaturedProblem(problem, Feature('linearly_transformed', options), 10, 3);
                kernel = optiprofiler_internal.FeatureKernel('linearly_transformed', options);
                [A, b, inv] = kernel.modifier_affine(3, problem);
                testCase.verifyTrue(isdiag(A) && cond(A) > 5, name);
                testCase.verifyFeatured(featured, problem, TestAffineStructure.expectedDiagonal(problem, A, b, inv), [name, ' unrotated']);

                options.rotated = true;
                featured = FeaturedProblem(problem, Feature('linearly_transformed', options), 10, 3);
                kernel = optiprofiler_internal.FeatureKernel('linearly_transformed', options);
                A = kernel.modifier_affine(3, problem);
                testCase.verifyFalse(isdiag(A), name);
                testCase.verifyFeatured(featured, problem, TestAffineStructure.expectedGeneric(problem, A, zeros(3, 1)), [name, ' rotated']);
            end
        end

        function strayEntryOfTheScalingsInverseKeepsBoundsAndOneOfItsMatrixIsPosed(testCase)
            % The framework builds both matrices itself, so the pair is
            % perturbed where it is produced (fixtures/affine); the code that
            % decides the representation is inherited unchanged.
            testCase.useAffineFixtures();
            options = struct('rotated', false, 'condition_factor', 4);
            corner = [0, 1, 0; 0, 0, 0; 0, 0, 0];
            % The matrix stays exactly diagonal, so the set is a box whatever
            % the inverse carries off its diagonal: at roundoff level, or above.
            perturbations = { ...
                @(A, inv) deal(A, inv + 0.25 * eps * min(abs(inv(1, 1)), abs(inv(2, 2))) * corner), ...
                @(A, inv) deal(A, inv + 1e-12 * abs(inv(1, 1)) * corner)};
            for t = 1:numel(perturbations)
                for k = 1:numel(TestAffineStructure.ProblemNames)
                    name = TestAffineStructure.ProblemNames{k};
                    label = sprintf('inverse perturbation %d %s', t, name);
                    problem = TestAffineStructure.makeProblem(name);
                    kernel = PerturbedAffineKernel('linearly_transformed', options, perturbations{t});
                    [A, b, inv] = kernel.modifier_affine(3, problem);
                    testCase.verifyTrue(isdiag(A) && ~isdiag(inv), label);  % the perturbation is in place
                    posed = TestAffineStructure.posedByKernel(kernel, 3, problem);
                    testCase.verifyStructure(posed, TestAffineStructure.expectedDiagonal(problem, A, b, inv), label);
                    testCase.verifyNoBoundIsLostOrDoubled(posed, problem, label);
                    testCase.verifyPosedIsScored(posed, A, b, @(y) problem.maxcv(A * y + b), problem, label);
                end
            end
            % A quarter of a unit in the last place of the diagonal of A is
            % still a coupling of two variables: it is posed, not ignored.
            coupling = @(A, inv) deal(A + 0.25 * eps * min(abs(A(1, 1)), abs(A(2, 2))) * corner, inv);
            for k = 1:numel(TestAffineStructure.ProblemNames)
                name = TestAffineStructure.ProblemNames{k};
                problem = TestAffineStructure.makeProblem(name);
                kernel = PerturbedAffineKernel('linearly_transformed', options, coupling);
                [A, b] = kernel.modifier_affine(3, problem);
                testCase.verifyFalse(isdiag(A), name);  % the perturbation is in place
                posed = TestAffineStructure.posedByKernel(kernel, 3, problem);
                testCase.verifyStructure(posed, TestAffineStructure.expectedGeneric(problem, A, b), ['coupling ', name]);
                testCase.verifyNoBoundIsLostOrDoubled(posed, problem, ['coupling ', name]);
                testCase.verifyPosedIsScored(posed, A, b, @(y) problem.maxcv(A * y + b), problem, ['coupling ', name]);
            end
        end

        function unusableConditionFactorsFailClosedAndOrdinaryOnesKeepWorking(testCase)
            problem = TestAffineStructure.makeProblem('linear');
            build = @(options) FeaturedProblem(problem, Feature('linearly_transformed', options), 10, 3);
            % 2^(+-1937) overflows: the matrix is not finite.
            testCase.verifyError(@() build(struct('rotated', false, 'condition_factor', 1e7)), 'MATLAB:Feature:AffineTransformationInvalid');
            testCase.verifyError(@() build(struct('rotated', true, 'condition_factor', 1e7)), 'MATLAB:Feature:AffineTransformationInvalid');
            % cond(A) = 2^95, unrotated. Every entry is exact and so is the
            % posed problem: a large condition number alone is no reason to
            % refuse.
            kernel = optiprofiler_internal.FeatureKernel('linearly_transformed', struct('rotated', false, 'condition_factor', 6000));
            [A, b, inv] = kernel.modifier_affine(3, problem);
            testCase.verifyGreaterThan(cond(A), 1e28);
            posed = TestAffineStructure.posedByKernel(kernel, 3, problem);
            testCase.verifyStructure(posed, TestAffineStructure.expectedDiagonal(problem, A, b, inv), 'cond 2^95');
            testCase.verifyNoBoundIsLostOrDoubled(posed, problem, 'cond 2^95');
            % The framework's own inverse is not held to the residual test of a
            % supplied one: a rotation with a large condition number is
            % accurate only to roundoff times that number, which is not an
            % inconsistency.
            for rotated = [false, true]
                for condition_factor = [0, 2, 5, 40]
                    featured = build(struct('rotated', rotated, 'condition_factor', condition_factor));
                    testCase.verifyFeatured(featured, problem, [], sprintf('rotated=%d condition_factor=%g', rotated, condition_factor));
                end
            end
            % cond(A) = 2^40 for n = 3. Both factors are built directly
            % (nothing is inverted), so each is accurate to roundoff, but their
            % product misses the identity by roundoff times the condition
            % number: more than a supplied inverse is allowed, and no
            % inconsistency. The posed rows need A only.
            kernel = optiprofiler_internal.FeatureKernel('linearly_transformed', ...
                struct('rotated', true, 'condition_factor', 2 * 40^2 / 3));
            [A, b, inv] = kernel.modifier_affine(3, problem);
            testCase.verifyGreaterThan(norm(A * inv - eye(3), 'fro'), 1e-8 * 3);
            posed = TestAffineStructure.posedByKernel(kernel, 3, problem);
            testCase.verifyStructure(posed, TestAffineStructure.expectedGeneric(problem, A, b), 'cond 2^40');
            testCase.verifyNoBoundIsLostOrDoubled(posed, problem, 'cond 2^40');
        end

        % ------------------------------------------------ composition, load and save

        function compositionsPoseTheScoredProblem(testCase)
            custom = @TestAffineStructure.customStage;
            scaling = struct('name', 'linearly_transformed', 'options', struct('rotated', false, 'condition_factor', 4));
            cases = { ...
                {custom(@TestAffineStructure.roundoffInverse), 'noisy'}, true; ...
                {'perturbed_x0', custom(@TestAffineStructure.roundoffMatrix), 'truncated'}, false; ...  % a coupling in A is posed
                {scaling, custom(@TestAffineStructure.roundoffInverse)}, true; ...
                {'permuted', custom(@TestAffineStructure.exactDiagonal)}, true; ...
                {custom(@TestAffineStructure.sloppyInverse), 'noisy'}, true; ...  % A is exactly diagonal: the set is a box
                {custom(@TestAffineStructure.dense), struct('name', 'linearly_transformed', 'options', struct('rotated', false))}, false};
            for c = 1:size(cases, 1)
                for k = 1:numel(TestAffineStructure.ProblemNames)
                    name = TestAffineStructure.ProblemNames{k};
                    label = sprintf('composition %d %s', c, name);
                    problem = TestAffineStructure.makeProblem(name);
                    featured = FeaturedProblem(problem, Feature(cases{c, 1}), 10, 3);
                    testCase.verifyEqual(featured.execution_strategy, 'composed-views', label);
                    testCase.verifyEqual(any(isfinite(featured.xl)) || any(isfinite(featured.xu)), cases{c, 2}, label);
                    testCase.verifyFeatured(featured, problem, [], label);
                    testCase.verifyEqual(featured.reference, TestAffineStructure.Reference, label);  % every stage here is a safe one
                end
            end
        end

        function savedTransformsKeepTheirStructureAndInvalidOnesStillFailClosed(testCase)
            problem = TestAffineStructure.makeProblem('linear');
            file = [tempname, '.mat'];
            testCase.addTeardown(@() deleteIfPresent(file));
            valid = {@TestAffineStructure.exactDiagonal, @TestAffineStructure.roundoffInverse, ...
                @TestAffineStructure.sloppyInverse, @TestAffineStructure.dense};
            for t = 1:numel(valid)
                label = func2str(valid{t});
                feature = Feature({TestAffineStructure.customStage(valid{t})});
                featured = FeaturedProblem(problem, Feature({TestAffineStructure.customStage(valid{t}), 'noisy'}), 10, 3);
                before = TestAffineStructure.posedOf(FeaturedProblem(problem, feature, 10, 3));
                save(file, 'feature', 'featured');
                loaded = load(file);
                % A saved specification poses the same problem after loading ...
                testCase.verifyEqual(TestAffineStructure.posedOf(FeaturedProblem(problem, loaded.feature, 10, 3)), before, label);
                % ... and a saved composed problem keeps its structure and its scoring.
                testCase.verifyEqual(TestAffineStructure.posedOf(loaded.featured), TestAffineStructure.posedOf(featured), label);
                testCase.verifyFeatured(loaded.featured, problem, [], label);
            end
            % A specification is data: it can be saved whatever its callbacks
            % return. The safeguard runs where the problem is built, so loading
            % cannot smuggle an invalid transform past it.
            invalid = { ...
                @TestAffineStructure.nanInMatrix, 'MATLAB:Feature:AffineTransformationInvalid'; ...
                @TestAffineStructure.wrongInverseShape, 'MATLAB:Feature:AffineTransformationNotInvertible'; ...
                @TestAffineStructure.materiallyWrongInverse, 'MATLAB:Feature:AffineTransformationNotInvertible'; ...
                @TestAffineStructure.numericallySingular, 'MATLAB:Feature:AffineTransformationNotInvertible'};
            for t = 1:size(invalid, 1)
                feature = Feature({TestAffineStructure.customStage(invalid{t, 1})});
                save(file, 'feature');
                loaded = load(file);
                testCase.verifyError(@() FeaturedProblem(problem, loaded.feature, 10, 3), invalid{t, 2}, func2str(invalid{t, 1}));
            end
        end

        % ------------------------------------------------ follow-up: the decision is exact

        function couplingBelowEveryToleranceIsPosedNotIgnored(testCase)
            % The base called 4e-16 negligible next to a diagonal of 1 (n * eps
            % is 4.4e-16) and posed the box [0, 1e16]^2 in the new variables.
            problem = Problem(struct('fun', @TestAffineStructure.sphere, 'x0', [1; 1], 'xl', [0; 0], ...
                'xu', [1e16; 1e16], 'reference', TestAffineStructure.Reference));
            [A, b] = TestAffineStructure.roundingLevelCoupling();
            stage = TestAffineStructure.customStage(@TestAffineStructure.roundingLevelCoupling);
            pipelines = {{stage}, {stage, 'noisy'}};
            % y = (-4, 1e16) is mapped to x = (0, 1e16), a vertex of the
            % original box, and that box in the new variables rejected it by 4;
            % it accepted (1e16, 1e16), which is mapped 4 outside.
            points = [-4, 1e16, 0, 3; 1e16, 1e16, 0, 5e15];
            violations = [0, 4, 0, 0];
            for c = 1:2
                label = sprintf('composed=%d', c - 1);
                featured = FeaturedProblem(problem, Feature(pipelines{c}), 10, 3);
                posed = TestAffineStructure.posedOf(featured);
                testCase.verifyStructure(posed, TestAffineStructure.expectedGeneric(problem, A, b), label);
                testCase.verifyNoBoundIsLostOrDoubled(posed, problem, label);
                for j = 1:numel(violations)
                    y = points(:, j);
                    testCase.verifyEqual(scoredViolation(featured, y), violations(j), 'AbsTol', 1e-9, label);
                    testCase.verifyEqual(TestAffineStructure.posedViolation(posed, y), violations(j), 'AbsTol', 1e-9, label);
                end
                testCase.verifyEqual(featured.reference, TestAffineStructure.Reference, label);
            end
        end

        function diagonalInverseThatIsNotTheReciprocalIsNotUsedForTheBounds(testCase)
            % The pair is accepted: it is consistent to 1e-9 and the contract
            % asks for 1e-8. The shortcut would scale the bounds by that
            % inverse; the rows of A are exact whatever the inverse is.
            [A, b, inv] = TestAffineStructure.inaccurateDiagonalInverse();
            stage = TestAffineStructure.customStage(@TestAffineStructure.inaccurateDiagonalInverse);
            pipelines = {{stage}, {stage, 'noisy'}};
            for c = 1:2
                for k = 1:numel(TestAffineStructure.ProblemNames)
                    name = TestAffineStructure.ProblemNames{k};
                    problem = TestAffineStructure.makeProblem(name);
                    featured = FeaturedProblem(problem, Feature(pipelines{c}), 10, 3);
                    testCase.verifyFeatured(featured, problem, TestAffineStructure.expectedGeneric(problem, A, b), ...
                        sprintf('composed=%d %s', c - 1, name));
                end
            end
            % Bounds of any size are posed to roundoff. The upper bounds that
            % diag(inv) gives: on the base the corner of the posed box, and
            % 1e7 outside the original one.
            problem = Problem(struct('fun', @TestAffineStructure.sphere, 'x0', [1; 1; 1], 'xl', zeros(3, 1), 'xu', 1e16 * ones(3, 1)));
            featured = FeaturedProblem(problem, Feature(pipelines{1}), 10, 3);
            scale = diag(inv);
            y = max(scale .* (problem.xl - b), scale .* (problem.xu - b));
            scored = scoredViolation(featured, y);
            testCase.verifyGreaterThan(scored, 1e6);
            testCase.verifyEqual(TestAffineStructure.posedViolation(TestAffineStructure.posedOf(featured), y), scored, 'RelTol', 1e-5);
        end

        % ------------------------------------------------ follow-up: one transform per problem

        function callbackWithAStateCannotMixTwoTransforms(testCase)
            % On the base every modifier asked the callback again. With the
            % dense answer for the bounds and the diagonal one for the rows,
            % the bounds became infinite and no bound row was added: every
            % finite bound was lost. The other way round every bound was posed
            % twice.
            testCase.useAffineFixtures();
            [diagonal_A, diagonal_b, diagonal_inv] = TestAffineStructure.exactDiagonal();
            for calls = 0:1
                for composed = [false, true]
                    for k = 1:numel(TestAffineStructure.ProblemNames)
                        name = TestAffineStructure.ProblemNames{k};
                        label = sprintf('calls=%d composed=%d %s', calls, composed, name);
                        problem = TestAffineStructure.makeProblem(name);
                        state = StatefulAffine('alternating', calls);
                        stages = TestAffineStructure.pipeline(@(s, p) state.transform(s, p), composed);
                        featured = FeaturedProblem(problem, Feature(stages), 10, 3);
                        testCase.verifyEqual(state.calls, calls + 1, [label, ' asked once']);
                        if calls == 0
                            expected = TestAffineStructure.expectedDiagonal(problem, diagonal_A, diagonal_b, diagonal_inv);
                        else
                            expected = TestAffineStructure.expectedGeneric(problem, TestAffineStructure.Dense, TestAffineStructure.B);
                        end
                        testCase.verifyFeatured(featured, problem, expected, label);
                        testCase.verifyEqual(state.calls, calls + 1, [label, ' and not again, 400 truth reads later']);
                        testCase.verifyEqual(featured.reference, TestAffineStructure.Reference, label);
                    end
                end
            end
        end

        function evaluationsUseTheTransformTheStructureWasBuiltWith(testCase)
            testCase.useAffineFixtures();
            problem = TestAffineStructure.makeProblem('linear');
            [A, b] = StatefulAffine.rotation(1);  % the first answer, and the only one that may be used
            points = -2 + 4 * rand(RandStream('mt19937ar', 'Seed', 1), 3, 6);
            for composed = [false, true]
                label = sprintf('composed=%d', composed);
                state = StatefulAffine('drifting');
                stages = TestAffineStructure.pipeline(@(s, p) state.transform(s, p), composed);
                featured = FeaturedProblem(problem, Feature(stages), 20, 3);
                testCase.verifyStructure(TestAffineStructure.posedOf(featured), TestAffineStructure.expectedGeneric(problem, A, b), label);
                testCase.verifyEqual(A * featured.x0 + b, problem.x0, 'AbsTol', 1e-14, label);
                testCase.verifyEqual(featured.fun_init, problem.fun(problem.x0), 'RelTol', 1e-13, label);
                for j = 1:size(points, 2)
                    y = points(:, j);
                    featured.fun(y);
                    testCase.verifyEqual(featured.fun_hist(end), problem.fun(A * y + b), 'RelTol', 1e-13, label);
                    testCase.verifyEqual(featured.maxcv_hist(end), problem.maxcv(A * y + b), 'AbsTol', 1e-13, label);
                    testCase.verifyEqual(featured.maxcv(y), problem.maxcv(A * y + b), 'AbsTol', 1e-13, label);
                end
                testCase.verifyEqual(state.calls, 1, label);
            end
        end

        function transformIsKeptPerSeedAndPerProblem(testCase)
            testCase.useAffineFixtures();
            first = TestAffineStructure.makeProblem('linear');
            second = TestAffineStructure.makeProblem('bounded');
            state = StatefulAffine('counting');
            kernel = optiprofiler_internal.FeatureKernel('custom', struct('mod_affine', @(s, p) state.transform(s, p)));
            for repeat = 1:3
                kernel.modifier_x0(3, first);
                TestAffineStructure.posedByKernel(kernel, 3, first);
                kernel.modifier_affine(3, first);
            end
            testCase.verifyEqual(state.calls, 1);
            kernel.modifier_affine(4, first);
            testCase.verifyEqual(state.calls, 2, 'another seed is another stream for the callback');
            kernel.modifier_bounds(4, second);
            testCase.verifyEqual(state.calls, 3, 'and another problem another question');
            kernel.modifier_linear_ub(4, second);
            testCase.verifyEqual(state.calls, 3);
            % A specification keeps no such state: every featured problem asks afresh.
            feature = Feature('custom', struct('mod_affine', @(s, p) state.transform(s, p)));
            FeaturedProblem(first, feature, 10, 3);
            FeaturedProblem(first, feature, 10, 3);
            testCase.verifyEqual(state.calls, 5);
        end

        function savedProblemKeepsTheTransformItWasBuiltWith(testCase)
            % A saved featured problem carries its structure. It has to carry
            % the map the structure was built with as well: asking the callback
            % again after loading may give another one.
            testCase.useAffineFixtures();
            problem = TestAffineStructure.makeProblem('linear');
            file = [tempname, '.mat'];
            testCase.addTeardown(@() deleteIfPresent(file));
            [A, b] = StatefulAffine.rotation(1);
            y = [0.3; -0.2; 0.7];
            for composed = [false, true]
                label = sprintf('composed=%d', composed);
                state = StatefulAffine('drifting');
                stages = TestAffineStructure.pipeline(@(s, p) state.transform(s, p), composed);
                featured = FeaturedProblem(problem, Feature(stages), 10, 3);
                save(file, 'featured');
                loaded = load(file);
                testCase.verifyStructure(TestAffineStructure.posedOf(loaded.featured), TestAffineStructure.expectedGeneric(problem, A, b), label);
                loaded.featured.fun(y);
                testCase.verifyEqual(loaded.featured.fun_hist(end), problem.fun(A * y + b), 'RelTol', 1e-13, label);
                testCase.verifyEqual(loaded.featured.maxcv(y), problem.maxcv(A * y + b), 'AbsTol', 1e-13, label);
            end
        end

        % ------------------------------------------------ follow-up: the inverse, from both sides

        function inverseHasToInvertFromBothSides(testCase)
            % Badly scaled, so that one product is the identity to 1e-16 while
            % the other misses it by 1: the base looked at the first only.
            problem = TestAffineStructure.makeProblem('bounded');
            not_invertible = 'MATLAB:Feature:AffineTransformationNotInvertible';
            transforms = {@TestAffineStructure.oneSidedInverse, @TestAffineStructure.oneSidedInverseMirror};
            for t = 1:numel(transforms)
                label = func2str(transforms{t});
                [A, ~, inv] = transforms{t}();
                testCase.verifyLessThanOrEqual(norm(A * inv - eye(3), 'fro'), 1e-8 * 3, label);
                testCase.verifyGreaterThan(norm(inv * A - eye(3), 'fro'), 0.9, label);
                kernel = optiprofiler_internal.FeatureKernel('custom', struct('mod_affine', transforms{t}));
                testCase.verifyError(@() kernel.modifier_affine(3, problem), not_invertible, label);
                % The transposed pairs fail on the other side, as they did before.
                kernel = optiprofiler_internal.FeatureKernel('custom', struct('mod_affine', @(s, p) deal(A', zeros(3, 1), inv')));
                testCase.verifyError(@() kernel.modifier_affine(3, problem), not_invertible, label);
            end
        end

        function neitherSetOfUnitsMakesAConsistentPairInconsistent(testCase)
            % A rotation with one variable in units of 1e13: new variable
            % (columns of A) or original one (rows). Roundoff times the ratio
            % of the units, 1e-3, is in one product or in the other; measured
            % against the terms the entries are summed from, both are roundoff.
            problem = TestAffineStructure.makeProblem('bounded');
            not_invertible = 'MATLAB:Feature:AffineTransformationNotInvertible';
            transforms = {@TestAffineStructure.columnScaledRotation, @TestAffineStructure.rowScaledRotation};
            large = [2, 1];
            for t = 1:numel(transforms)
                label = func2str(transforms{t});
                [A, b, inv] = transforms{t}();
                residuals = [norm(A * inv - eye(3), 'fro'), norm(inv * A - eye(3), 'fro')];
                testCase.verifyGreaterThan(residuals(large(t)), 1e-8 * 3, label);  % the plain rule refuses one side
                testCase.verifyLessThan(residuals(3 - large(t)), 1e-12, label);
                kernel = optiprofiler_internal.FeatureKernel('custom', struct('mod_affine', transforms{t}));
                [kept_A, kept_b, kept_inv] = kernel.modifier_affine(3, problem);
                testCase.verifyEqual({kept_A, kept_b, kept_inv}, {A, b, inv}, label);
            end
            % What was inconsistent stays inconsistent: the scale of the terms
            % never excuses an error of the size of the terms.
            transforms = {@TestAffineStructure.singular, @TestAffineStructure.identityAsInverse, @TestAffineStructure.materiallyWrongInverse};
            for t = 1:numel(transforms)
                [A, b, inv] = transforms{t}();
                pairs = {A, inv; A', inv'};
                for side = 1:2
                    kernel = optiprofiler_internal.FeatureKernel('custom', struct('mod_affine', @(s, p) deal(pairs{side, 1}, b, pairs{side, 2})));
                    testCase.verifyError(@() kernel.modifier_affine(3, problem), not_invertible, func2str(transforms{t}));
                end
            end
            [A, b, inv] = TestAffineStructure.rowScaledRotation();
            inv(1, 2) = inv(1, 2) * (1 + 1e-6);  % one entry of the inverse off by 1e-6 of its size
            kernel = optiprofiler_internal.FeatureKernel('custom', struct('mod_affine', @(s, p) deal(A, b, inv)));
            testCase.verifyError(@() kernel.modifier_affine(3, problem), not_invertible);
        end

        % ------------------------------------------------ follow-up: nothing finite leaves the range

        function nothingFiniteLeavesTheRangeInTheTransport(testCase)
            % xu - b with 1e308 + 1e308, bub - aub * b with aub * b = 2e318,
            % aub * A with 1e200 * 1e200, inv * (x0 - b). The scale guard of
            % the base looked at the scaled product only, and an infinite
            % factor is "no bound" to it. The base posed infinite upper bounds,
            % rows with an infinite right-hand side (which count as no
            % constraint), right-hand sides of -Inf and NaN, and coefficients
            % of Inf.
            bounds = 'MATLAB:Feature:AffineBoundsNotRepresentable';
            rows = 'MATLAB:Feature:AffineLinearConstraintsNotRepresentable';
            start = 'MATLAB:Feature:AffineInitialPointNotRepresentable';
            make = @(varargin) Problem(struct('fun', @TestAffineStructure.sphere, 'x0', zeros(3, 1), varargin{:}));
            far_box = make('xl', zeros(3, 1), 'xu', 1e308 * ones(3, 1));
            cases = { ...
                far_box, @TestAffineStructure.shiftedOutOfRange, bounds; ...
                far_box, @TestAffineStructure.denseShiftedOutOfRange, bounds; ...
                make('xl', [1e308; 0; 0], 'xu', [1e308; 1; 1]), @TestAffineStructure.denseShiftedOutOfRange, bounds; ...  % a fixed variable
                make('aub', [1e10, 1e10, 0], 'bub', 1), @TestAffineStructure.shiftedTheOtherWay, rows; ...
                make('aeq', [1e10, -3e10, 0], 'beq', 2), @TestAffineStructure.shiftedTheOtherWay, rows; ...
                make('aub', [1e200, 1, 0], 'bub', 1), @TestAffineStructure.hugeScaling, rows; ...
                make('aeq', [1e200, 0, 1], 'beq', 0), @TestAffineStructure.hugeScaling, rows; ...
                Problem(struct('fun', @TestAffineStructure.sphere, 'x0', 1e308 * ones(3, 1))), @TestAffineStructure.shiftedOutOfRange, start};
            for c = 1:size(cases, 1)
                stage = TestAffineStructure.customStage(cases{c, 2});
                label = sprintf('case %d %s', c, func2str(cases{c, 2}));
                testCase.verifyError(@() FeaturedProblem(cases{c, 1}, Feature({stage}), 10, 3), cases{c, 3}, label);
                testCase.verifyError(@() FeaturedProblem(cases{c, 1}, Feature({stage, 'noisy'}), 10, 3), cases{c, 3}, label);
            end
            % With supplied inequality rows the inequality modifier of the
            % framework does not run (and nothing of its is replaced here: the
            % only finite bounds are those of the fixed variable). The equality
            % of that variable is still posed by the framework, and has its own
            % guard.
            only_fixed = make('xl', [1e308; -Inf; -Inf], 'xu', [1e308; Inf; Inf]);
            stage = struct('name', 'custom', 'options', struct('mod_affine', @TestAffineStructure.denseShiftedOutOfRange, ...
                'mod_linear_ub', @(s, p) deal([1, 1, 0], 5)));
            testCase.verifyError(@() FeaturedProblem(only_fixed, Feature({stage}), 10, 3), bounds, 'fixed variable, supplied rows');
            testCase.verifyError(@() FeaturedProblem(only_fixed, Feature({stage, 'noisy'}), 10, 3), bounds, 'fixed variable, supplied rows');
            % The scaling feature is guarded as well: diag(A) reaches 2^47 for
            % this condition factor.
            feature = Feature('linearly_transformed', struct('rotated', false, 'condition_factor', 6000));
            testCase.verifyError(@() FeaturedProblem(make('aub', [1e300, 0, 1e300], 'bub', 1), feature, 10, 3), rows);
            testCase.verifyError(@() FeaturedProblem(Problem(struct('fun', @TestAffineStructure.sphere, 'x0', [1e300; 0; 0])), feature, 10, 3), start);
        end

        function boundThatIsInfiniteAlreadyStaysAsItIs(testCase)
            % An infinite bound is no bound: nothing finite is lost, so nothing
            % is refused, however large the finite data next to it.
            problem = Problem(struct('fun', @TestAffineStructure.sphere, 'x0', [0.5; 0.5; 0.5], ...
                'xl', [-Inf; -1e300; -Inf], 'xu', [Inf; Inf; 1e300], 'aub', [1, 1, 0], 'bub', 1e300));
            transforms = {@TestAffineStructure.exactDiagonal, @TestAffineStructure.dense};
            for t = 1:numel(transforms)
                [A, b, inv] = transforms{t}();
                if t == 1
                    expected = TestAffineStructure.expectedDiagonal(problem, A, b, inv);
                else
                    expected = TestAffineStructure.expectedGeneric(problem, A, b);
                end
                stage = TestAffineStructure.customStage(transforms{t});
                pipelines = {{stage}, {stage, 'noisy'}};
                for c = 1:2
                    label = sprintf('%s composed=%d', func2str(transforms{t}), c - 1);
                    featured = FeaturedProblem(problem, Feature(pipelines{c}), 10, 3);
                    testCase.verifyStructure(TestAffineStructure.posedOf(featured), expected, label);
                    testCase.verifyTrue(all(isfinite(featured.bub)) && all(isfinite(featured.aub(:))), label);
                end
            end
        end

        % ------------------------------------------------ follow-up: the class of the data

        function dataOfAProblemAreDoublesWhateverClassTheyWereGivenIn(testCase)
            % MATLAB combines an integer with a SCALAR double in integer
            % arithmetic, rounded and saturated, and refuses to combine it with
            % an array. With int32 data the base posed the box [-1, 2] for the
            % exact [-0.25, 1.75] and started at 1 for 0.25; and the truth, in
            % the same arithmetic, called x = 4.3 feasible for xu = 4.
            sphere = @(x) sum(double(x(:)).^2);
            problem = Problem(struct('fun', sphere, 'x0', int32(1), 'xl', int32(0), 'xu', int32(4)));
            testCase.verifyClass(problem.x0, 'double');
            testCase.verifyClass(problem.xl, 'double');
            testCase.verifyClass(problem.xu, 'double');
            testCase.verifyEqual(problem.maxcv(4.3), 0.3, 'AbsTol', 1e-12);
            featured = FeaturedProblem(problem, Feature('custom', struct('mod_affine', @(s, p) deal(2, 0.5, 0.5))), 10, 3);
            testCase.verifyEqual([featured.xl, featured.xu, featured.x0], [-0.25, 1.75, 0.25]);
            testCase.verifyEqual(featured.maxcv(1.9), 0.3, 'AbsTol', 1e-12);  % x = 2 * 1.9 + 0.5 = 4.3

            % Several variables, several classes: MATLAB:mixedClasses on the base.
            integers = struct('fun', sphere, 'x0', int32([1; 1; 1]), 'xl', int8([-1; -2; -3]), 'xu', uint16([3; 5; 4]), ...
                'aub', int32([1, 1, 0]), 'bub', int32(5), 'aeq', single([1, 0, -1]), 'beq', single(0.25));
            doubles = struct('fun', sphere, 'x0', [1; 1; 1], 'xl', [-1; -2; -3], 'xu', [3; 5; 4], ...
                'aub', [1, 1, 0], 'bub', 5, 'aeq', [1, 0, -1], 'beq', 0.25);
            names = {'x0', 'xl', 'xu', 'aub', 'bub', 'aeq', 'beq'};
            given = Problem(integers);
            for k = 1:numel(names)
                testCase.verifyClass(given.(names{k}), 'double', names{k});
                testCase.verifyEqual(given.(names{k}), double(integers.(names{k})), names{k});
            end
            transforms = {@TestAffineStructure.exactDiagonal, @TestAffineStructure.dense};
            for t = 1:numel(transforms)
                stage = TestAffineStructure.customStage(transforms{t});
                pipelines = {{stage}, {stage, 'noisy'}};
                for c = 1:2
                    label = sprintf('%s composed=%d', func2str(transforms{t}), c - 1);
                    actual = FeaturedProblem(Problem(integers), Feature(pipelines{c}), 10, 3);
                    expected = FeaturedProblem(Problem(doubles), Feature(pipelines{c}), 10, 3);
                    testCase.verifyEqual(TestAffineStructure.posedOf(actual), TestAffineStructure.posedOf(expected), label);
                    testCase.verifyEqual(actual.x0, expected.x0, label);
                end
            end

            % Single precision data are transported in double precision, from
            % the numbers that were stored: 2e-8 of the bound on the base.
            problem = Problem(struct('fun', sphere, 'x0', single(1), 'xl', single(0.1), 'xu', single(4.7)));
            featured = FeaturedProblem(problem, Feature('custom', struct('mod_affine', @(s, p) deal(3, 0.1, 1 / 3))), 10, 3);
            testCase.verifyClass(featured.xu, 'double');
            testCase.verifyEqual(featured.xu, (1 / 3) * (double(single(4.7)) - 0.1));
        end

        % ------------------------------------------------ second follow-up: nothing finite is lost to underflow

        function nothingFiniteIsLostToUnderflow(testCase)
            % The finding: xl = 1e-200, xu = 2e-200, A = 1e200, inv = 1e-200.
            % Both scaled bounds are 1e-400 and 2e-400, which round to 0: the
            % base posed the single point y = 0, which is mapped to x = 0,
            % outside the interval. The same loss removed a coefficient row
            % (aub * A with 2^-600 * 2^-600: a row of zeros is no constraint),
            % a right-hand side (0 - 2^-1200) and the initial point (the base
            % started at 0, where log(x1) is -Inf, for a finite fun(x0)).
            bounds = 'MATLAB:Feature:AffineBoundsNotRepresentable';
            rows = 'MATLAB:Feature:AffineLinearConstraintsNotRepresentable';
            start = 'MATLAB:Feature:AffineInitialPointNotRepresentable';
            sphere = @TestAffineStructure.sphere;
            make = @(varargin) Problem(struct('fun', sphere, 'x0', zeros(3, 1), varargin{:}));
            huge = @(s, p) deal(1e200, 0, 1e-200);
            shift = @(s, p) deal(eye(3), [2^-600; 0; 0], eye(3));
            logarithm = Problem(struct('fun', @(x) log(x(1)) + x(2)^2 + x(3)^2, 'x0', [2^-600; 1; 1]));
            testCase.verifyTrue(isfinite(logarithm.fun(logarithm.x0)));
            cases = { ...
                Problem(struct('fun', sphere, 'x0', 0, 'xl', 1e-200, 'xu', 2e-200)), huge, bounds; ...         % the finding
                Problem(struct('fun', sphere, 'x0', 1.5e-200, 'xl', 1e-200, 'xu', 2e-200)), huge, start; ...   % from inside the interval the point is lost first
                make('xl', [2^-500; -1; -1], 'xu', [1; 1; 1]), TestAffineStructure.powerScaling(1023 - 500), bounds; ...   % half of realmin
                make('xl', [2^-500; -1; -1], 'xu', [1; 1; 1]), TestAffineStructure.powerScaling(1050 - 500), bounds; ...   % a subnormal number
                make('xl', [2^-500; -1; -1], 'xu', [1; 1; 1]), TestAffineStructure.powerScaling(1100 - 500), bounds; ...   % zero
                make('aub', [2^-600, 0, 0], 'bub', 1), TestAffineStructure.powerScaling(-600), rows; ...
                make('aeq', [2^-600, 0, 0], 'beq', 0), TestAffineStructure.powerScaling(-600), rows; ...
                make('aub', [2^-600, 0, 0], 'bub', 0), shift, rows; ...
                make('aeq', [2^-600, 0, 0], 'beq', 0), shift, rows; ...
                logarithm, TestAffineStructure.powerScaling(600), start};
            for c = 1:size(cases, 1)
                for composed = [false, true]
                    label = sprintf('case %d composed=%d', c, composed);
                    feature = Feature(TestAffineStructure.pipeline(cases{c, 2}, composed));
                    testCase.verifyError(@() FeaturedProblem(cases{c, 1}, feature, 10, 3), cases{c, 3}, label);
                end
            end
            try
                FeaturedProblem(cases{1, 1}, Feature(TestAffineStructure.pipeline(huge, false)), 10, 3);
            catch exception
                testCase.verifySubstring(exception.message, 'underflows');
            end
            % The scaling feature is guarded as well: diag(inv) reaches 2^-47
            % for this condition factor, and 1e-300 becomes a subnormal number.
            feature = Feature('linearly_transformed', struct('rotated', false, 'condition_factor', 6000));
            testCase.verifyError(@() FeaturedProblem(make('xl', [-1; -1; 1e-300], 'xu', [1; 1; 1]), feature, 10, 3), bounds);
            testCase.verifyError(@() FeaturedProblem(make('aub', [1e-300, 0, 0], 'bub', 1), feature, 10, 3), rows);
        end

        function underflowBoundaryAndWhatIsNotUnderflow(testCase)
            sphere = @TestAffineStructure.sphere;
            for composed = [false, true]
                label = sprintf('composed=%d', composed);
                % The boundary is the smallest normal number: 2^-500 scaled by
                % 2^-522 is exactly realmin, and it is posed.
                problem = Problem(struct('fun', sphere, 'x0', zeros(3, 1), 'xl', [2^-500; -1; -1], 'xu', [1; 1; 1]));
                featured = FeaturedProblem(problem, Feature(TestAffineStructure.pipeline(TestAffineStructure.powerScaling(1022 - 500), composed)), 10, 3);
                testCase.verifyEqual(featured.xl(1), realmin, label);
                % 0 times anything is 0 exactly: nothing is lost, nothing refused.
                problem = Problem(struct('fun', sphere, 'x0', zeros(3, 1), 'xl', [0; -1; -1], 'xu', [Inf; 1; 1]));
                featured = FeaturedProblem(problem, Feature(TestAffineStructure.pipeline(TestAffineStructure.powerScaling(900), composed)), 10, 3);
                testCase.verifyEqual(featured.xl, [0; -1; -1], label);
                % Next to a right-hand side that is a normal number, a shift
                % that underflows is below its rounding.
                shift = @(s, p) deal(eye(3), [2^-600; 0; 0], eye(3));
                problem = Problem(struct('fun', sphere, 'x0', zeros(3, 1), 'aub', [2^-600, 0, 0], 'bub', 5));
                featured = FeaturedProblem(problem, Feature(TestAffineStructure.pipeline(shift, composed)), 10, 3);
                testCase.verifyEqual(featured.bub, 5, label);
                % Cancellation is not underflow: an entry that is zero because
                % its terms cancel has lost nothing, its terms are normal.
                A = [2, 1, 0; 1, 1, 0; 0, 0, 1];
                mixing = @(s, p) deal(A, [1; 1; 0], A \ eye(3));
                problem = Problem(struct('fun', sphere, 'x0', [1; 1; 0], 'aub', [1, -1, 0], 'bub', 0, 'aeq', [1, -2, 0], 'beq', -1));
                featured = FeaturedProblem(problem, Feature(TestAffineStructure.pipeline(mixing, composed)), 10, 3);
                testCase.verifyEqual(featured.aub, [1, 0, 0], label);
                testCase.verifyEqual(featured.bub, 0, label);
                testCase.verifyEqual(featured.aeq, [0, -1, 0], label);
                testCase.verifyEqual(featured.beq, 0, label);
                testCase.verifyEqual(featured.x0, zeros(3, 1), label);
            end
        end

        % ------------------------------------------------ second follow-up: derivatives in the solver's variables

        function gradientOfTheFinding(testCase)
            % f(x) = ||x||^2, A = diag(2, 1), y = (1, 1): the function the
            % solver sees is 4 * y1^2 + y2^2 with gradient (8, 2). The base
            % returned grad(y) = (2, 2), the gradient of another function at
            % another point.
            problem = Problem(struct('fun', @TestAffineStructure.sphere, 'x0', [1; 1], 'grad', @(x) 2 * x(:), 'hess', @(x) 2 * eye(2)));
            scaling = @(s, p) deal(diag([2, 1]), [0; 0], diag([0.5, 1]));
            featured = FeaturedProblem(problem, Feature('custom', struct('mod_affine', scaling)), 100, 3);
            testCase.verifyEqual(featured.grad([1; 1]), [8; 2]);
            testCase.verifyEqual(featured.hess([1; 1]), [8, 0; 0, 2]);
            testCase.verifyEqual(featured.grad([1; 1]), TestAffineStructure.centralDifferences(@(y) featured.fun(y), [1; 1])', 'RelTol', 1e-6);
        end

        function everyDerivativeFollowsTheChainRule(testCase)
            problem = TestAffineStructure.smoothProblem(true);
            features = TestAffineStructure.variableChanges();
            Y = [0.2; -0.3; 0.4];
            fd = @TestAffineStructure.centralDifferences;
            close = {'RelTol', 1e-6, 'AbsTol', 1e-6};
            exact = {'RelTol', 1e-13, 'AbsTol', 1e-13};
            for k = 1:numel(features)
                label = features{k}.name;
                if k > 1 && strcmp(label, features{k - 1}.name), label = [label, ' (diagonal)']; end %#ok<AGROW>
                seed = 5;
                [A, b] = TestAffineStructure.transformationOf(features{k}, seed, problem);
                while isequal(A, eye(3)), seed = seed + 1; [A, b] = TestAffineStructure.transformationOf(features{k}, seed, problem); end
                featured = FeaturedProblem(problem, features{k}, 10000, seed);
                x = A * Y + b;
                % Exactly the chain rule of x = A * y + b ...
                testCase.verifyEqual(featured.grad(Y), A' * problem.grad(x), exact{:}, label);
                testCase.verifyEqual(featured.hess(Y), A' * problem.hess(x) * A, exact{:}, label);
                testCase.verifyEqual(featured.jcub(Y), problem.jcub(x) * A, exact{:}, label);
                testCase.verifyEqual(featured.jceq(Y), problem.jceq(x) * A, exact{:}, label);
                returned = {featured.hcub(Y), featured.hceq(Y)}; original = {problem.hcub(x), problem.hceq(x)};
                for c = 1:2
                    testCase.verifyEqual(numel(returned{c}), numel(original{c}), label);
                    for i = 1:numel(original{c})
                        testCase.verifyEqual(returned{c}{i}, A' * original{c}{i} * A, exact{:}, label);
                    end
                end
                % ... which is what the functions the solver evaluates vary by.
                testCase.verifyEqual(featured.grad(Y), fd(@(y) featured.fun(y), Y)', close{:}, label);
                testCase.verifyEqual(featured.hess(Y), fd(@(y) featured.grad(y), Y), close{:}, label);
                testCase.verifyEqual(featured.jcub(Y), fd(@(y) featured.cub(y), Y), close{:}, label);
                testCase.verifyEqual(featured.jceq(Y), fd(@(y) featured.ceq(y), Y), close{:}, label);
                H = featured.hcub(Y);
                for i = 1:numel(H)
                    testCase.verifyEqual(H{i}, fd(@(y) TestAffineStructure.rowOf(featured.jcub(y), i), Y), close{:}, label);
                end
                H = featured.hceq(Y);
                testCase.verifyEqual(H{1}, fd(@(y) featured.jceq(y), Y), close{:}, label);
            end
        end

        function derivativesCostNothingStayAbsentAndCheckThePoint(testCase)
            testCase.useAffineFixtures();
            Y = [0.2; -0.3; 0.4];
            names = {'grad', 'hess', 'jcub', 'jceq', 'hcub', 'hceq'};
            % They cost nothing and read the kept transformation.
            state = StatefulAffine('counting');
            featured = FeaturedProblem(TestAffineStructure.smoothProblem(true), Feature('custom', struct('mod_affine', @(s, p) state.transform(s, p))), 10, 3);
            for k = 1:numel(names), featured.(names{k})(Y); end
            testCase.verifyEqual(state.calls, 1);
            testCase.verifyEqual([featured.n_eval_fun, featured.n_eval_cub, featured.n_eval_ceq], [0, 0, 0]);
            testCase.verifyEmpty(featured.fun_hist);
            testCase.verifyEmpty(featured.maxcv_hist);
            features = TestAffineStructure.variableChanges();
            without = TestAffineStructure.smoothProblem(false);
            with = TestAffineStructure.smoothProblem(true);
            for f = 1:numel(features)
                % An absent derivative stays absent.
                featured = FeaturedProblem(without, features{f}, 10, 5);
                for k = 1:numel(names)
                    testCase.verifyEmpty(featured.(names{k})(Y), sprintf('%s %s', features{f}.name, names{k}));
                end
                % The point is checked as before, by Problem.
                featured = FeaturedProblem(with, features{f}, 10, 5);
                for k = 1:numel(names)
                    testCase.verifyError(@() featured.(names{k})([0.1; 0.2]), ['MATLAB:Problem:WrongSizeInputFor', upper(names{k})], names{k});
                end
            end
        end

        function withoutAChangeOfVariablesDerivativesAreThoseOfTheProblemAndACompositionHasNone(testCase)
            % The established passthrough: these features change values or the
            % initial point, never the variables, and the derivative methods
            % describe the callbacks of the problem, not the observed values.
            problem = TestAffineStructure.smoothProblem(true);
            Y = [0.2; -0.3; 0.4];
            names = {'grad', 'hess', 'jcub', 'jceq', 'hcub', 'hceq'};
            features = {'plain', 'noisy', 'truncated', 'perturbed_x0', 'random_nan', 'quantized', ...
                'unrelaxable_constraints', 'nonquantifiable_constraints'};
            for f = 1:numel(features)
                featured = FeaturedProblem(problem, Feature(features{f}), 10, 5);
                for k = 1:numel(names)
                    testCase.verifyEqual(featured.(names{k})(Y), problem.(names{k})(Y), sprintf('%s %s', features{f}, names{k}));
                end
            end
            featured = FeaturedProblem(problem, Feature({TestAffineStructure.customStage(@TestAffineStructure.denseWithShift), 'noisy'}), 10, 3);
            for k = 1:numel(names)
                testCase.verifyError(@() featured.(names{k})(Y), 'MATLAB:FeaturedProblem:UnsupportedCompositeDerivative', names{k});
            end
        end

        % ------------------------------------------------ second follow-up: the initial point

        function pointIsVerifiedWhereItIsUsed(testCase)
            % A = I with inv(1, 2) = 1e-12 is an identity to 1e-12 from both
            % sides, far inside every tolerance on the matrices, and it moves
            % the feasible x0 = (0, 1e14) to (100, 1e14): the base started 99
            % outside the bounds. No such tolerance bounds an error at a point.
            sphere = @TestAffineStructure.sphere;
            stray = @(s, p) deal(eye(2), [0; 0], [1, 1e-12; 0, 1]);
            problem = Problem(struct('fun', sphere, 'x0', [0; 1e14], 'xl', [-1; 0], 'xu', [1; 1e14]));
            for composed = [false, true]
                featured = FeaturedProblem(problem, Feature(TestAffineStructure.pipeline(stray, composed)), 10, 3);
                testCase.verifyEqual(featured.x0, [0; 1e14]);
                testCase.verifyEqual(featured.maxcv_init, 0);
                testCase.verifyEqual(featured.fun_init, sphere(problem.x0));
            end
            % The test is by component: next to a component of 1e14 an error of
            % 1e-3 in the other one is 1e-17 of the norm, and it is not excused.
            kernel = optiprofiler_internal.FeatureKernel('custom', struct('mod_affine', @(s, p) deal(eye(2), [0; 0], [1, 1e-17; 0, 1])));
            start = Problem(struct('fun', sphere, 'x0', [1; 1e14]));
            testCase.verifyEqual(kernel.modifier_x0(3, start), [1; 1e14]);
            % Within the allowance the supplied inverse is believed: one unit in the last place.
            nudged = [1 + eps, 0; 0, 1];
            kernel = optiprofiler_internal.FeatureKernel('custom', struct('mod_affine', @(s, p) deal(eye(2), [0; 0], nudged)));
            testCase.verifyEqual(kernel.modifier_x0(3, start), nudged * [1; 1e14]);
        end

        function goodInverseGivesThePointItAlwaysGaveAndEveryPointIsMappedBack(testCase)
            % Bitwise: the framework's own inverse, and a supplied one that is
            % an inverse to roundoff at the point.
            problem = TestAffineStructure.makeProblem('linear');
            options = {struct('rotated', false, 'condition_factor', 4), struct('rotated', true, 'condition_factor', 40), ...
                struct('rotated', true, 'condition_factor', 6000)};
            for k = 1:numel(options)
                feature = Feature('linearly_transformed', options{k});
                [~, ~, inv] = TestAffineStructure.transformationOf(feature, 3, problem);
                featured = FeaturedProblem(problem, feature, 10, 3);
                testCase.verifyEqual(featured.x0, inv * problem.x0, sprintf('scaling %d', k));
            end
            transforms = {@TestAffineStructure.exactDiagonal, @TestAffineStructure.dense, @TestAffineStructure.sloppyInverse, ...
                @TestAffineStructure.badlyScaledRotation, @TestAffineStructure.rowScaledRotation, @TestAffineStructure.columnScaledRotation};
            for t = 1:numel(transforms)
                label = func2str(transforms{t});
                [A, b, inv] = transforms{t}();
                for k = 1:numel(TestAffineStructure.ProblemNames)
                    problem = TestAffineStructure.makeProblem(TestAffineStructure.ProblemNames{k});
                    featured = FeaturedProblem(problem, Feature('custom', struct('mod_affine', transforms{t})), 10, 3);
                    allowed = 64 * 3 * eps * (abs(A) * abs(featured.x0) + abs(b) + abs(problem.x0));
                    testCase.verifyTrue(all(abs(A * featured.x0 + b - problem.x0) <= allowed), label);
                    if ~isequal(transforms{t}, @TestAffineStructure.sloppyInverse)
                        testCase.verifyEqual(featured.x0, inv * (problem.x0 - b), label);
                    end
                end
            end
        end

        function allowanceIsTheRoundingOfTheEvaluationAtThatPoint(testCase)
            % x = (y1 + 1e15 * y2, y2): at x0 = (1/3, 1/7) the first component
            % is evaluated as 1.4e14 - 1.4e14 + 1/3, which double precision
            % resolves to 0.03, whatever y is. The exact inverse is therefore
            % accepted and its point is mapped back within the allowance, which
            % is NOT the rounding of x0: it is what the truth, which goes
            % through the same map at every point, can resolve there.
            A = [1, 1e15; 0, 1]; inv = [1, -1e15; 0, 1]; x0 = [1 / 3; 1 / 7];
            kernel = optiprofiler_internal.FeatureKernel('custom', struct('mod_affine', @(s, p) deal(A, [0; 0], inv)));
            point = kernel.modifier_x0(3, Problem(struct('fun', @TestAffineStructure.sphere, 'x0', x0)));
            testCase.verifyEqual(point, inv * x0);
            error_ = abs(A * point - x0);
            testCase.verifyTrue(all(error_ <= 64 * 2 * eps * (abs(A) * abs(point) + abs(x0))));
            testCase.verifyLessThanOrEqual(error_(1), 0.0625);
        end

        % ------------------------------------------------ second follow-up: decisions made explicit

        function consistencyAllowsTheRoundingOfTheProductsAndNothingMore(testCase)
            % 330805f measured each residual against max(1, terms), which
            % allowed 1e-8 OF THE TERMS. For this matrix the terms are 2e15, and
            % an inverse with one entry off by 1e-9 of its size was accepted:
            % the initial point was pulled back to a point that is mapped 1e6
            % away from x0. 3a43a19 refused it. What the terms justify is the
            % rounding of the products.
            not_invertible = 'MATLAB:Feature:AffineTransformationNotInvertible';
            problem = Problem(struct('fun', @TestAffineStructure.sphere, 'x0', [1; 1]));
            A = [1, 1e15; 0, 1]; exact = [1, -1e15; 0, 1]; sloppy = [1, -1e15 + 1e6; 0, 1];
            featured = FeaturedProblem(problem, Feature('custom', struct('mod_affine', @(s, p) deal(A, [0; 0], exact))), 10, 3);
            testCase.verifyEqual(A * featured.x0, problem.x0);
            testCase.verifyError(@() FeaturedProblem(problem, Feature('custom', struct('mod_affine', @(s, p) deal(A, [0; 0], sloppy))), 10, 3), not_invertible);
            testCase.verifyError(@() FeaturedProblem(problem, Feature('custom', struct('mod_affine', @(s, p) deal(A', [0; 0], sloppy'))), 10, 3), not_invertible);
            % A pair that is consistent to roundoff is accepted at any condition
            % number: built as linearly_transformed builds its own, D * Q' and
            % Q * D^-1, with a condition number of 1e12, each product misses the
            % identity by 1e-5, a few units of rounding of its terms.
            [Q, ~] = qr(randn(RandStream('mt19937ar', 'Seed', 7), 3));
            d = [1e-6; 1; 1e6];
            matrix = diag(d) * Q'; inverse = Q * diag(1 ./ d);
            testCase.verifyGreaterThan(max(norm(matrix * inverse - eye(3), 'fro'), norm(inverse * matrix - eye(3), 'fro')), 1e-8 * 3);
            kernel = optiprofiler_internal.FeatureKernel('custom', struct('mod_affine', @(s, p) deal(matrix, zeros(3, 1), inverse)));
            kept = kernel.modifier_affine(3, TestAffineStructure.makeProblem('bounded'));
            testCase.verifyEqual(kept, matrix);
        end

        function modBoundsReplacesTheBoundsAndNothingElse(testCase)
            % A supplied modifier replaces its own component, verbatim. The
            % bounds of the problem live in the bounds if the map keeps them
            % there, and then mod_bounds replaces them; under any other map they
            % are linear rows of the framework, which mod_bounds does not touch.
            % Never fewer constraints than before; the reference is unknown
            % either way.
            box = @(s, p) deal(-9 * ones(3, 1), 9 * ones(3, 1));
            problem = TestAffineStructure.makeProblem('linear');
            for composed = [false, true]
                label = sprintf('composed=%d', composed);
                tail = {}; if composed, tail = {'noisy'}; end
                build = @(transform) FeaturedProblem(problem, Feature([{struct('name', 'custom', 'options', ...
                    struct('mod_affine', transform, 'mod_bounds', box))}, tail]), 10, 3);
                featured = build(@TestAffineStructure.exactDiagonal);
                testCase.verifyEqual(featured.xl, -9 * ones(3, 1), label);
                testCase.verifyEqual(featured.xu, 9 * ones(3, 1), label);
                testCase.verifyEqual(featured.aub, problem.aub * diag(TestAffineStructure.D), label);  % no row of a bound
                testCase.verifyEmpty(featured.reference, label);
                featured = build(@TestAffineStructure.dense);
                expected = TestAffineStructure.expectedGeneric(problem, TestAffineStructure.Dense, TestAffineStructure.B);
                testCase.verifyEqual(featured.xl, -9 * ones(3, 1), label);
                testCase.verifyEqual(featured.xu, 9 * ones(3, 1), label);
                testCase.verifyEqual(featured.aub, expected.aub, label);  % the bounds of the problem, as rows
                testCase.verifyEqual(featured.bub, expected.bub, label);
                testCase.verifyEmpty(featured.reference, label);
            end
        end

        function integerBeyond2To53IsRoundedAndDeprecatedConveniencesKeepNothing(testCase)
            % Conversion to double precision is rounding to nearest, for an
            % integer as for a decimal literal, the same in both paths and both
            % languages. The rounded triple is the one that is validated, kept
            % and used, by the structure and by the truth alike.
            testCase.useAffineFixtures();
            big = int64(2)^53 + 1;
            A = int64(eye(3)); A(1, 1) = big;
            inverse = diag([1 / double(big), 1, 1]);
            problem = TestAffineStructure.makeProblem('bounded');
            y = [5e-16; 0.5; 0.5];  % mapped to x1 = 4.5, which violates xu = 3 by 1.5
            for composed = [false, true]
                label = sprintf('composed=%d', composed);
                integers = FeaturedProblem(problem, Feature(TestAffineStructure.pipeline(@(s, p) deal(A, zeros(3, 1, 'int64'), inverse), composed)), 10, 3);
                rounded = FeaturedProblem(problem, Feature(TestAffineStructure.pipeline(@(s, p) deal(diag([2^53, 1, 1]), zeros(3, 1), inverse), composed)), 10, 3);
                testCase.verifyEqual(TestAffineStructure.posedOf(integers), TestAffineStructure.posedOf(rounded), label);
                testCase.verifyEqual(integers.x0, rounded.x0, label);
                testCase.verifyEqual(integers.maxcv(y), problem.maxcv([2^53 * 5e-16; 0.5; 0.5]), label);
                testCase.verifyGreaterThan(integers.maxcv(y), 1, label);
            end
            % Each call of a deprecated Feature.modifier_* builds a kernel of
            % its own, so a callback with a state is asked again by each of
            % them: one transformation per problem is a property of
            % FeaturedProblem.
            state = StatefulAffine('counting');
            feature = Feature('custom', struct('mod_affine', @(s, p) state.transform(s, p)));
            testCase.verifyWarning(@() feature.modifier_bounds(3, problem), 'MATLAB:Feature:DeprecatedModifier');
            testCase.verifyWarning(@() feature.modifier_linear_ub(3, problem), 'MATLAB:Feature:DeprecatedModifier');
            testCase.verifyEqual(state.calls, 2);
        end
    end
end


function cv = scoredViolation(featured, y)
    % The violation the truth channel records for the solver point y.
    [~, cv] = featured.evaluateTruth(y);
end


function deleteIfPresent(file)
    if exist(file, 'file') == 2
        delete(file);
    end
end
