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

        function stage = customStage(transform)
            stage = struct('name', 'custom', 'options', struct('mod_affine', transform));
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

        function roundoffInEitherMatrixKeepsBoundsAsBounds(testCase)
            % The defect. With roundoff in the inverse the bounds became
            % infinite and no bound row was added: the bounds were gone. With
            % roundoff in the matrix the bounds stayed AND were added again as rows.
            transforms = {@TestAffineStructure.roundoffInverse, @TestAffineStructure.roundoffMatrix};
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

        function pairNotDiagonalToRoundoffTakesTheGenericPath(testCase)
            % Not roundoff, not inconsistent: the two matrices do not tell the
            % same structural story, so the representation that needs only A
            % is used. On the base a coupling in the inverse lost every bound,
            % and one in the matrix posed every bound twice.
            transforms = {@TestAffineStructure.sloppyInverse, @TestAffineStructure.sloppyMatrix, ...
                @TestAffineStructure.couplingHiddenByScale};
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
                @TestAffineStructure.numericallySingular, not_invertible, not_invertible, 'numerically singular'};
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

        function roundoffInEitherMatrixOfTheScalingKeepsBoundsAsBounds(testCase)
            % The framework builds both matrices itself, so the pair is
            % perturbed where it is produced (fixtures/affine); the code that
            % decides the representation is inherited unchanged.
            old_path = path;
            testCase.addTeardown(@() path(old_path));
            addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'fixtures', 'affine'));
            options = struct('rotated', false, 'condition_factor', 4);
            perturbations = { ...
                @(A, inv) deal(A, inv + 0.25 * eps * min(abs(inv(1, 1)), abs(inv(2, 2))) * [0, 1, 0; 0, 0, 0; 0, 0, 0]), ...
                @(A, inv) deal(A + 0.25 * eps * min(abs(A(1, 1)), abs(A(2, 2))) * [0, 1, 0; 0, 0, 0; 0, 0, 0], inv)};
            for t = 1:numel(perturbations)
                for k = 1:numel(TestAffineStructure.ProblemNames)
                    name = TestAffineStructure.ProblemNames{k};
                    label = sprintf('perturbation %d %s', t, name);
                    problem = TestAffineStructure.makeProblem(name);
                    kernel = PerturbedAffineKernel('linearly_transformed', options, perturbations{t});
                    [A, b, inv] = kernel.modifier_affine(3, problem);
                    testCase.verifyFalse(isdiag(A) && isdiag(inv), label);  % the perturbation is in place
                    posed = TestAffineStructure.posedByKernel(kernel, 3, problem);
                    testCase.verifyStructure(posed, TestAffineStructure.expectedDiagonal(problem, A, b, inv), label);
                    testCase.verifyNoBoundIsLostOrDoubled(posed, problem, label);
                    testCase.verifyPosedIsScored(posed, A, b, @(y) problem.maxcv(A * y + b), problem, label);
                end
            end
            % An inverse that is not diagonal to roundoff: the generic path, and no bound lost.
            sloppy = @(A, inv) deal(A, inv + 1e-12 * abs(inv(1, 1)) * [0, 1, 0; 0, 0, 0; 0, 0, 0]);
            for k = 1:numel(TestAffineStructure.ProblemNames)
                name = TestAffineStructure.ProblemNames{k};
                problem = TestAffineStructure.makeProblem(name);
                kernel = PerturbedAffineKernel('linearly_transformed', options, sloppy);
                [A, b] = kernel.modifier_affine(3, problem);
                posed = TestAffineStructure.posedByKernel(kernel, 3, problem);
                testCase.verifyStructure(posed, TestAffineStructure.expectedGeneric(problem, A, b), ['sloppy ', name]);
                testCase.verifyNoBoundIsLostOrDoubled(posed, problem, ['sloppy ', name]);
                testCase.verifyPosedIsScored(posed, A, b, @(y) problem.maxcv(A * y + b), problem, ['sloppy ', name]);
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
                {'perturbed_x0', custom(@TestAffineStructure.roundoffMatrix), 'truncated'}, true; ...
                {scaling, custom(@TestAffineStructure.roundoffInverse)}, true; ...
                {'permuted', custom(@TestAffineStructure.exactDiagonal)}, true; ...
                {custom(@TestAffineStructure.sloppyInverse), 'noisy'}, false; ...
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
