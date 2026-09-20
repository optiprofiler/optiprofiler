classdef FeatureKernel < handle
%FEATUREKERNEL Internal per-runtime adapter for the established MATLAB kernels.
% Input is an already-normalized atomic stage; no specification or experiment
% validation/defaulting occurs here. The numerical methods are extracted from
% ac4a67e Feature unchanged, including legacy read order and RNG mixing.
%
% PAYLOAD_MIXER selects how the per-query random stream is seeded from the
% run/stage seed and the observed payload (values, point, served index):
%   'legacy-product' (default): the established default_rng mixer of the
%       identity and single-feature strategies (seed_policy legacy-run-seed).
%       It multiplies the payload, so any zero element removes the payload.
%   'horner32-words': the composed-views mixer of seed_policy
%       matlab-stage-horner32-v2 (see horner32_payload_rng).
% Construction streams (initial point, permutation, rotation, custom
% structure) take no payload and use default_rng under both settings.
    properties (SetAccess = private)
        name
        options
        payload_mixer = 'legacy-product'
    end
    properties (Access = private)
        % The affine transformation of this runtime, produced once: a struct
        % with the fields problem, seed, A, b and inv (see modifier_affine). It
        % is saved with the kernel, so that a loaded featured problem goes on
        % with the map its structure was built with.
        kept_affine = []
    end
    methods
        function obj = FeatureKernel(name, options, payload_mixer)
            obj.name = name;
            obj.options = options;
            if nargin > 2
                if ~ismember(payload_mixer, {'legacy-product', 'horner32-words'})
                    error('MATLAB:Feature:UnknownPayloadMixer', 'Unknown payload mixer: %s.', payload_mixer);
                end
                obj.payload_mixer = payload_mixer;
            end
        end
        function x0 = modifier_x0(obj, seed, problem)
            %{
            Modify the initial point.

            Parameters
            ----------
            seed : int
                Seed used to generate random numbers.
            problem : Problem
                Problem for which the initial point is modified.

            Returns
            -------
            x0 : double, size (n,)
                Modified initial point.
            %}

            switch obj.name
                case FeatureName.CUSTOM.value
                    % If the user specifies a custom modifier for the initial point, use it.
                    if isfield(obj.options, FeatureOptionKey.MOD_X0.value)
                        rand_stream_custom = obj.default_rng(seed);
                        x0 = obj.options.(FeatureOptionKey.MOD_X0.value)(rand_stream_custom, problem);
                        return;
                    end
                    % If the user does not specify a custom modifier for the initial point but specifies
                    % a custom affine transformation, we need to apply the inverse of the affine
                    % transformation to the initial point.
                    if isfield(obj.options, FeatureOptionKey.MOD_AFFINE.value)
                        [A, b, inv] = obj.modifier_affine(seed, problem);
                        x0 = pulledBack(A, inv, problem.x0, b);
                    else
                        x0 = problem.x0;
                    end
                case FeatureName.PERTURBED_X0.value
                    % Use max(1, norm(x0)) to avoid no perturbation when x0 is zero.
                    rand_stream_perturbed_x0 = obj.default_rng(seed);
                    perturbation_level = obj.options.(FeatureOptionKey.PERTURBATION_LEVEL.value) * max(1, norm(problem.x0));
                    if strcmp(obj.options.(FeatureOptionKey.DISTRIBUTION.value), 'gaussian')
                        x0 = problem.x0 + perturbation_level * rand_stream_perturbed_x0.randn(problem.n, 1);
                    elseif strcmp(obj.options.(FeatureOptionKey.DISTRIBUTION.value), 'spherical')
                        perturbation = rand_stream_perturbed_x0.randn(problem.n, 1);
                        x0 = problem.x0 + perturbation_level * perturbation / norm(perturbation);
                    else
                        x0 = problem.x0 + perturbation_level * obj.options.(FeatureOptionKey.DISTRIBUTION.value)(rand_stream_perturbed_x0, problem.n);
                    end
                case FeatureName.PERMUTED.value
                    % Note that we need to apply the reverse permutation to the initial point so that
                    % the new problem is mathematically equivalent to the original one.
                    rand_stream_permuted = obj.default_rng(seed);
                    permutation = rand_stream_permuted.randperm(problem.n);
                    [~, reverse_permutation] = sort(permutation);
                    x0 = problem.x0(reverse_permutation);
                case FeatureName.LINEARLY_TRANSFORMED.value
                    % Apply the inverse of the affine transformation to the initial point.
                    [A, ~, inv] = obj.modifier_affine(seed, problem);
                    x0 = pulledBack(A, inv, problem.x0, []);
                otherwise
                    x0 = problem.x0;
            end
        end

        function tf = changesVariables(obj)
            % Whether the stage is a change of variables x = A * y + b other
            % than the identity by construction.
            tf = ismember(obj.name, {FeatureName.PERMUTED.value, FeatureName.LINEARLY_TRANSFORMED.value}) || ...
                (strcmp(obj.name, FeatureName.CUSTOM.value) && isfield(obj.options, FeatureOptionKey.MOD_AFFINE.value));
        end

        function [A, b, inv] = modifier_affine(obj, seed, problem)
            %{
            Generate an invertible matrix A and a vector b for the affine transformation applied to the variables.

            Parameters
            ----------
            seed : int
                Seed used to generate random numbers.
            problem : Problem
                Problem for which the affine transformation is generated.
            Returns
            -------
            A : double, size (n, n)
                Matrix of the affine transformation.
            b : double, size (n, 1)
                Vector of the affine transformation.
            inv : double, size (n, n)
                Inverse of the matrix A.

            The triple is produced once per problem and seed and kept: the
            initial point, the bounds, both kinds of linear constraints and
            every evaluation read the same one. A mod_affine with a state (a
            counter, the global random stream) may answer differently when
            asked again, and the bounds of one map with the linear constraints
            of another are the structure of no problem. A specification keeps
            nothing: a kernel serves one trial.
            %}

            kept = obj.kept_affine;
            if ~isempty(kept) && kept.problem == problem && isequal(kept.seed, seed)
                A = kept.A;
                b = kept.b;
                inv = kept.inv;
                return;
            end
            [A, b, inv] = obj.generateAffine(seed, problem);
            obj.kept_affine = struct('problem', problem, 'seed', seed, 'A', A, 'b', b, 'inv', inv);
        end
    end

    methods (Access = private)
        function [A, b, inv] = generateAffine(obj, seed, problem)
            % The triple of modifier_affine, produced anew: the only place
            % that calls mod_affine.

            % Default values
            A = eye(problem.n);
            b = zeros(problem.n, 1);
            inv = eye(problem.n);

            switch obj.name
                case FeatureName.CUSTOM.value
                    if isfield(obj.options, FeatureOptionKey.MOD_AFFINE.value)
                        rand_stream_custom = obj.default_rng(seed);
                        [A, b, inv] = obj.options.(FeatureOptionKey.MOD_AFFINE.value)(rand_stream_custom, problem);
                        % The inverse is supplied by user code: validate the
                        % triple, including that A * inv is an identity matrix.
                        [A, b, inv] = checkedAffine(A, b, inv, problem.n, true);
                    end
                case FeatureName.PERMUTED.value
                    % Generate a random permutation matrix.
                    rand_stream_permuted = obj.default_rng(seed);
                    permutation = rand_stream_permuted.randperm(problem.n);
                    A = eye(problem.n);
                    A = A(permutation, :);
                    inv = A';
                    b = zeros(problem.n, 1);
                case FeatureName.LINEARLY_TRANSFORMED.value
                    % Generate A in the form D * Q' with D being a diagonal matrix and Q being an
                    % orthogonal matrix.

                    % Generate a random rotation matrix Q if the 'rotated' option is set to true.
                    if obj.options.(FeatureOptionKey.ROTATED.value)
                        %{
                        We generate a random orthogonal matrix Q following the uniform distribution
                        on O(n). The method refers to a note written by the late Professor Nicholas
                        Higham, see:
                        https://nhigham.com/2020/04/22/what-is-a-random-orthogonal-matrix/
                        and an answer from MathStackExchange, see:
                        https://math.stackexchange.com/a/4891933/1088047
                        %}
                        rand_stream_linearly_transformed = obj.default_rng(seed);
                        [Q, R] = qr(rand_stream_linearly_transformed.randn(problem.n));
                        Q(:, diag(R) < 0) = -Q(:, diag(R) < 0);
                    else
                        Q = eye(problem.n);
                    end
                    % The extreme exponents differ by sqrt(condition_factor*n/2),
                    % so cond(A)=2^sqrt(condition_factor*n/2) for n>=2. Orthogonal
                    % rotation preserves singular values; n=1 has cond(A)=1.
                    log_condition_number = sqrt(obj.options.(FeatureOptionKey.CONDITION_FACTOR.value) * problem.n / 2);
                    power = linspace(-log_condition_number/2, log_condition_number/2, problem.n);
                    A = diag(2.^power) * Q';
                    inv = Q * diag(2.^-power);
                    % Built here, so consistent by construction; a huge
                    % condition factor can still overflow or be singular to
                    % working precision.
                    [A, b, inv] = checkedAffine(A, b, inv, problem.n, false);
                otherwise
                    % Do nothing
            end
        end
    end

    methods
        function [xl, xu] = modifier_bounds(obj, seed, problem)
            %{
            Modify the bounds.

            Parameters
            ----------
            seed : int
                Seed used to generate random numbers.
            problem : Problem
                Problem for which the bounds are modified.

            Returns
            -------
            xl : double, size (n,)
                Modified lower bounds.
            xu : double, size (n,)
                Modified upper bounds.
            %}

            switch obj.name
                case FeatureName.CUSTOM.value
                    % If the user specifies a custom modifier for the bounds, use it.
                    if isfield(obj.options, FeatureOptionKey.MOD_BOUNDS.value)
                        rand_stream_custom = obj.default_rng(seed);
                        [xl, xu] = obj.options.(FeatureOptionKey.MOD_BOUNDS.value)(rand_stream_custom, problem);
                        return;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.MOD_AFFINE.value)
                        xl = problem.xl;
                        xu = problem.xu;
                        return;
                    end
                    % If the user does not specify a custom modifier for the bounds but specifies a
                    % custom affine transformation, we need to specially handle the bounds.
                    [A, b, inv] = obj.modifier_affine(seed, problem);
                    if ~affineIsDiagonal(A, inv)
                        % Generic representation: the bounds are posed as
                        % linear rows (see modifier_linear_ub and
                        % modifier_linear_eq, which read the same decision).
                        xl = -Inf(problem.n, 1);
                        xu = Inf(problem.n, 1);
                        return;
                    end
                    % Diagonal shortcut: the bounds stay bounds, scaled by the inverse.
                    [xl, xu] = scaledBounds(diag(inv), shifted(problem.xl, b), shifted(problem.xu, b));
                case FeatureName.PERMUTED.value
                    % Note that we need to apply the reverse permutation to the bounds so that the new
                    % problem is mathematically equivalent to the original one.
                    rand_stream_permuted = obj.default_rng(seed);
                    permutation = rand_stream_permuted.randperm(problem.n);
                    [~, reverse_permutation] = sort(permutation);
                    xl = problem.xl(reverse_permutation);
                    xu = problem.xu(reverse_permutation);
                case FeatureName.LINEARLY_TRANSFORMED.value
                    % Apply the inverse of the affine transformation to the bounds.
                    [A, ~, inv] = obj.modifier_affine(seed, problem);
                    if ~affineIsDiagonal(A, inv)
                        xl = -Inf(problem.n, 1);
                        xu = Inf(problem.n, 1);
                        return;
                    end
                    [xl, xu] = scaledBounds(diag(inv), problem.xl, problem.xu);
                otherwise
                    xl = problem.xl;
                    xu = problem.xu;
            end
        end

        function [aub, bub] = modifier_linear_ub(obj, seed, problem)
            %{
            Modify the linear inequality constraints.

            Parameters
            ----------
            seed : int
                Seed used to generate random numbers.
            problem : Problem
                Problem for which the linear inequality constraints are modified.

            Returns
            -------
            aub : double, size (m_linear_ub, n)
                Modified matrix of linear inequality constraints.
            bub : double, size (m_linear_ub,)
                Modified right-hand side vector of linear inequality
                constraints.
            %}

            switch obj.name
                case FeatureName.CUSTOM.value
                    % If the user specifies a custom modifier for the linear inequality constraints,
                    % use it.
                    if isfield(obj.options, FeatureOptionKey.MOD_LINEAR_UB.value)
                        refuseToReplaceBoundRows(obj, seed, problem, 'mod_linear_ub', false);
                        rand_stream_custom = obj.default_rng(seed);
                        [aub, bub] = obj.options.(FeatureOptionKey.MOD_LINEAR_UB.value)(rand_stream_custom, problem);
                        return;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.MOD_AFFINE.value)
                        aub = problem.aub;
                        bub = problem.bub;
                        return;
                    end
                    % If the user does not specify a custom modifier for the linear inequality
                    % constraints but specifies a custom affine transformation, we need to specially
                    % handle the linear inequality constraints.
                    [A, b, inv] = obj.modifier_affine(seed, problem);
                    [aub, bub] = composedRows(problem.aub, problem.bub, A, b);
                    if affineIsDiagonal(A, inv)  % the bounds stayed bounds (modifier_bounds read the same decision)
                        return;
                    end
                    %{
                    We need to specially handle bound constraints and linear inequality constraints.

                    Bound constraints
                    xl <= A * x + b <= xu
                    should be modified to
                    A * x <= xu - b
                    and
                    -A * x <= -xl + b
                    when A is not diagonal.

                    Linear inequality constraints
                    turn out to be
                    (aub * A) * x <= bub - aub * b
                    %}

                    % Pick out the indices of lower bounds who are not -Inf and upper bounds who are
                    % not Inf since later we will not transform them into linear inequality constraints.
                    idx_lb = ~isinf(problem.xl);
                    idx_ub = ~isinf(problem.xu);
                    % Remove the indices, of which the lower and upper bound are equal, since later
                    % we will put them into the linear equality constraints.
                    idx_eq = find(problem.xl == problem.xu);
                    idx_lb(idx_eq) = false;
                    idx_ub(idx_eq) = false;
                    upper = shifted(problem.xu, b);
                    lower = shifted(problem.xl, b);
                    if isempty(problem.aub)
                        aub = [A(idx_ub, :); -A(idx_lb, :)];
                        bub = [upper(idx_ub); -lower(idx_lb)];
                        return;
                    end
                    aub = [A(idx_ub, :); -A(idx_lb, :); aub];
                    bub = [upper(idx_ub); -lower(idx_lb); bub];
                case FeatureName.PERMUTED.value
                    rand_stream_permuted = obj.default_rng(seed);
                    permutation = rand_stream_permuted.randperm(problem.n);
                    [~, reverse_permutation] = sort(permutation);
                    aub = problem.aub(:, reverse_permutation);
                    bub = problem.bub;
                case FeatureName.LINEARLY_TRANSFORMED.value
                    % Similar to the case in the custom feature where a custom affine transformation
                    % is specified.
                    [A, ~, inv] = obj.modifier_affine(seed, problem);
                    [aub, bub] = composedRows(problem.aub, problem.bub, A, []);
                    if affineIsDiagonal(A, inv)
                        return;
                    end
                    idx_lb = ~isinf(problem.xl);
                    idx_ub = ~isinf(problem.xu);
                    idx_eq = find(problem.xl == problem.xu);
                    idx_lb(idx_eq) = false;
                    idx_ub(idx_eq) = false;
                    if isempty(problem.aub)
                        aub = [A(idx_ub, :); -A(idx_lb, :)];
                        bub = [problem.xu(idx_ub); -problem.xl(idx_lb)];
                        return;
                    end
                    aub = [A(idx_ub, :); -A(idx_lb, :); aub];
                    bub = [problem.xu(idx_ub); -problem.xl(idx_lb); bub];
                otherwise
                    aub = problem.aub;
                    bub = problem.bub;
            end
        end

        function [aeq, beq] = modifier_linear_eq(obj, seed, problem)
            %{
            Modify the linear equality constraints.

            Parameters
            ----------
            seed : int
                Seed used to generate random numbers.
            problem : Problem
                Problem for which the linear equality constraints are modified.

            Returns
            -------
            aeq : double, size (m_linear_eq, n)
                Modified matrix of linear equality constraints.
            beq : double, size (m_linear_eq,)
                Modified right-hand side vector of linear equality constraints.
            %}

            switch obj.name
                case FeatureName.CUSTOM.value
                    % If the user specifies a custom modifier for the linear equality constraints, use it.
                    if isfield(obj.options, FeatureOptionKey.MOD_LINEAR_EQ.value)
                        refuseToReplaceBoundRows(obj, seed, problem, 'mod_linear_eq', true);
                        rand_stream_custom = obj.default_rng(seed);
                        [aeq, beq] = obj.options.(FeatureOptionKey.MOD_LINEAR_EQ.value)(rand_stream_custom, problem);
                        return;
                    end
                    if ~isfield(obj.options, FeatureOptionKey.MOD_AFFINE.value)
                        aeq = problem.aeq;
                        beq = problem.beq;
                        return;
                    end
                    % If the user does not specify a custom modifier for the linear equality constraints but specifies a custom affine transformation, we need to specially handle the linear equality constraints.
                    [A, b, inv] = obj.modifier_affine(seed, problem);
                    [aeq, beq] = composedRows(problem.aeq, problem.beq, A, b);
                    if affineIsDiagonal(A, inv)  % the bounds stayed bounds (modifier_bounds read the same decision)
                        return;
                    end
                    %{
                    We need to specially handle bound constraints and linear equality constraints.

                    Bound constraints
                    xl <= A * x + b <= xu
                    with xl = xu should be modified to
                    A * x = xu - b
                    when A is not diagonal.

                    Linear equality constraints
                    turn out to be
                    (aeq * A) * x = beq - aeq * b
                    %}

                    % Pick out the indices, of which the lower and upper bound are equal.
                    idx_eq = find(problem.xl == problem.xu);
                    fixed = shifted(problem.xu, b);
                    if isempty(problem.aeq)
                        aeq = A(idx_eq, :);
                        beq = fixed(idx_eq);
                        return;
                    end
                    aeq = [A(idx_eq, :); aeq];
                    beq = [fixed(idx_eq); beq];
                case FeatureName.PERMUTED.value
                    rand_stream_permuted = obj.default_rng(seed);
                    permutation = rand_stream_permuted.randperm(problem.n);
                    [~, reverse_permutation] = sort(permutation);
                    aeq = problem.aeq(:, reverse_permutation);
                    beq = problem.beq;
                case FeatureName.LINEARLY_TRANSFORMED.value
                    % Similar to the case in the custom feature where a custom affine transformation is specified.
                    [A, ~, inv] = obj.modifier_affine(seed, problem);
                    [aeq, beq] = composedRows(problem.aeq, problem.beq, A, []);
                    if affineIsDiagonal(A, inv)
                        return;
                    end
                    idx_eq = find(problem.xl == problem.xu);
                    if isempty(problem.aeq)
                        aeq = A(idx_eq, :);
                        beq = problem.xu(idx_eq);
                        return;
                    end
                    aeq = [A(idx_eq, :); aeq];
                    beq = [problem.xu(idx_eq); beq];
                otherwise
                    aeq = problem.aeq;
                    beq = problem.beq;
            end
        end

        function f = modifier_fun(obj, x, seed, problem, n_eval)
            %{
            Modify the objective function value.

            Parameters
            ----------
            x : double, size (n,)
                Point at which the objective function is evaluated.
            seed : int
                Seed used to generate random numbers.
            problem : Problem
                Optimization problem to be modified.
            n_eval : int
                Number of evaluations of the objective function.
                (We will use it to generate random streams so that evaluating
                the same point multiple times will not lead to the same
                random numbers.)

            Returns
            -------
            f : double
                Modified objective function value.
            %}

            % Convert x into a cell array. We will later use it to generate
            % random streams so that randomness of each point is independent.
            xCell = num2cell(x);

            f = problem.fun(x);

            switch obj.name
                case FeatureName.CUSTOM.value
                    if isfield(obj.options, FeatureOptionKey.MOD_FUN.value)
                        rand_stream_custom = obj.payloadStream(seed, f, xCell{:}, n_eval);
                        f = obj.options.(FeatureOptionKey.MOD_FUN.value)(x, rand_stream_custom, problem);
                        return;
                    end
                case FeatureName.NOISY.value
                    noise = obj.computeNoise(x, seed, n_eval, f, 1);
                    f = obj.applyNoise(f, noise);
                case FeatureName.RANDOM_NAN.value
                    rand_stream_random_nan = obj.payloadStream(seed, f, xCell{:}, n_eval);
                    if rand_stream_random_nan.rand() < obj.options.(FeatureOptionKey.NAN_RATE.value)
                        f = NaN;
                    end
                case FeatureName.TRUNCATED.value
                    if isnan(f) || isinf(f)
                        % If f is NaN or Inf, we do not need to truncate it.
                        % Note that if f is NaN or Inf, digits will be set to NaN or Inf respectively, which will lead
                        % to an error when calling 'round(f, digits)'.
                        return;
                    end
                    rand_stream_truncated = obj.payloadStream(seed, f, xCell{:}, n_eval);
                    if f == 0
                        digits = obj.options.(FeatureOptionKey.SIGNIFICANT_DIGITS.value) - 1;
                    else
                        digits = obj.options.(FeatureOptionKey.SIGNIFICANT_DIGITS.value) - floor(log10(abs(f))) - 1;
                    end
                    f = round(f, digits);
                    % Round f to the desired number of significant digits. (We can also use 'fix(x)'
                    % to round towards zero.) We can shorten the above code by using 'round(f, digits, "significant")'
                    % directly, but we want to keep the same format as the Python code. The code before
                    % 'round(f, digits)' cannot be removed since the default choice of 'round' in MATLAB
                    % will consider digits in relation to the decimal point.
                    if obj.options.(FeatureOptionKey.PERTURBED_TRAILING_DIGITS.value)
                        if f >= 0
                            f = f + rand_stream_truncated.rand() * 10 ^ (-digits);
                        else
                            f = f - rand_stream_truncated.rand() * 10 ^ (-digits);
                        end
                    end
                case FeatureName.UNRELAXABLE_CONSTRAINTS.value
                    [~, maxcv_bounds, maxcv_linear, maxcv_nonlinear] = problem.maxcv(x, true);
                    if obj.options.(FeatureOptionKey.UNRELAXABLE_BOUNDS.value) && maxcv_bounds > 0
                        f = Inf;
                    elseif obj.options.(FeatureOptionKey.UNRELAXABLE_LINEAR_CONSTRAINTS.value) && maxcv_linear > 0
                        f = Inf;
                    elseif obj.options.(FeatureOptionKey.UNRELAXABLE_NONLINEAR_CONSTRAINTS.value) && maxcv_nonlinear > 0
                        f = Inf;
                    end
                case FeatureName.QUANTIZED.value
                    mesh_size = obj.options.(FeatureOptionKey.MESH_SIZE.value);
                    if strcmp(obj.options.(FeatureOptionKey.MESH_TYPE.value), 'relative')
                        mesh_size = mesh_size .* max(1, abs(x));
                    end
                    x = mesh_size .* round(x ./ mesh_size);  % round(x) rounds x to the closest integer.
                    f = problem.fun(x);
                otherwise
                    % Do nothing
            end
        end

        function cub_ = modifier_cub(obj, x, seed, problem, n_eval_cub)
            %{
            Modify the values of the nonlinear inequality constraints.

            Parameters
            ----------
            x : double, size (n,)
                Point at which the nonlinear inequality constraints are evaluated.
            seed : int
                Seed used to generate random numbers.
            problem : Problem
                Problem for which the nonlinear inequality constraints are modified.
            n_eval_cub : int
                Number of evaluations of the nonlinear inequality constraints.
                (We will use it to generate random streams so that evaluating
                the same point multiple times will not lead to the same
                random numbers.)

            Returns
            -------
            cub_ : double, size (m_nonlinear_ub,)
                Modified values of the nonlinear inequality constraints.
            %}

            % Convert x into a cell array. We will later use it to generate
            % random streams so that randomness of each point is independent.
            xCell = num2cell(x);

            cub_ = problem.cub(x);
            % If cub_ is empty, return directly!
            if isempty(cub_)
                return;
            end
            cubCell = num2cell(cub_);

            switch obj.name
                case FeatureName.CUSTOM.value
                    if isfield(obj.options, FeatureOptionKey.MOD_CUB.value)
                        rand_stream_custom = obj.payloadStream(seed, cubCell{:}, xCell{:}, n_eval_cub);
                        cub_ = obj.options.(FeatureOptionKey.MOD_CUB.value)(x, rand_stream_custom, problem);
                        return;
                    end
                case FeatureName.NOISY.value
                    % Similar to the case in the modifier_fun method.
                    noise = obj.computeNoise(x, seed, n_eval_cub, cubCell, size(cub_));
                    cub_ = obj.applyNoise(cub_, noise);
                case FeatureName.RANDOM_NAN.value
                    % Similar to the case in the modifier_fun method.
                    rand_stream_random_nan = obj.payloadStream(seed, cubCell{:}, xCell{:}, n_eval_cub);
                    cub_(rand_stream_random_nan.rand(size(cub_)) < obj.options.(FeatureOptionKey.NAN_RATE.value)) = NaN;
                case FeatureName.TRUNCATED.value
                    % Similar to the case in the modifier_fun method.
                    rand_stream_truncated = obj.payloadStream(seed, cubCell{:}, xCell{:}, n_eval_cub);
                    digits = zeros(size(cub_));
                    digits(cub_ == 0) = obj.options.(FeatureOptionKey.SIGNIFICANT_DIGITS.value) - 1;
                    digits(cub_ ~= 0) = obj.options.(FeatureOptionKey.SIGNIFICANT_DIGITS.value) - floor(log10(abs(cub_(cub_ ~= 0)))) - 1;
                    % Note that `floor(log10(abs(NaN))) = NaN` and `floor(log10(abs(Inf))) = Inf` in MATLAB.
                    for i_cub = 1:numel(cub_)
                        if ~isnan(cub_(i_cub)) && ~isinf(cub_(i_cub))
                            cub_(i_cub) = round(cub_(i_cub), digits(i_cub));
                        end
                    end
                    if obj.options.(FeatureOptionKey.PERTURBED_TRAILING_DIGITS.value)
                        cub_(cub_ >= 0) = cub_(cub_ >= 0) + rand_stream_truncated.rand(size(cub_(cub_ >= 0))) .* (10 .^ (-digits(cub_ >= 0)));
                        cub_(cub_ < 0) = cub_(cub_ < 0) - rand_stream_truncated.rand(size(cub_(cub_ < 0))) .* (10 .^ (-digits(cub_ < 0)));
                    end
                case FeatureName.NONQUANTIFIABLE_CONSTRAINTS.value
                    % Set the elements whose value are less than or equal to 0 to 0.
                    cub_(cub_ <= 0) = 0;
                    % Preserve NaN: an undefined oracle is not a known violation.
                    cub_(cub_ > 0) = 1;
                case FeatureName.QUANTIZED.value
                    % Similar to the case in the modifier_fun method.
                    mesh_size = obj.options.(FeatureOptionKey.MESH_SIZE.value);
                    if strcmp(obj.options.(FeatureOptionKey.MESH_TYPE.value), 'relative')
                        mesh_size = mesh_size .* max(1, abs(x));
                    end
                    x = mesh_size .* round(x ./ mesh_size);
                    cub_ = problem.cub(x);
                otherwise
                    % Do nothing
            end
        end

        function ceq_ = modifier_ceq(obj, x, seed, problem, n_eval_ceq)
            %{
            Modify the values of the nonlinear equality constraints.

            Parameters
            ----------
            x : double, size (n,)
                Point at which the nonlinear equality constraints are evaluated.
            seed : int
                Seed used to generate random numbers.
            problem : Problem
                Problem for which the nonlinear equality constraints are modified.
            n_eval_ceq : int
                Number of evaluations of the nonlinear equality constraints.
                (We will use it to generate random streams so that evaluating
                the same point multiple times will not lead to the same
                random numbers.)

            Returns
            -------
            ceq_ : double, size (m_nonlinear_eq,)
                Modified values of the nonlinear equality constraints.
            %}

            % Convert x into a cell array. We will later use it to generate
            % random streams so that randomness of each point is independent.
            xCell = num2cell(x);

            ceq_ = problem.ceq(x);
            % If ceq_ is empty, return directly!
            if isempty(ceq_)
                return;
            end
            ceqCell = num2cell(ceq_);

            switch obj.name
                case FeatureName.CUSTOM.value
                    if isfield(obj.options, FeatureOptionKey.MOD_CEQ.value)
                        rand_stream_custom = obj.payloadStream(seed, ceqCell{:}, xCell{:}, n_eval_ceq);
                        ceq_ = obj.options.(FeatureOptionKey.MOD_CEQ.value)(x, rand_stream_custom, problem);
                        return;
                    end
                case FeatureName.NOISY.value
                    % Similar to the case in the modifier_fun method.
                    noise = obj.computeNoise(x, seed, n_eval_ceq, ceqCell, size(ceq_));
                    ceq_ = obj.applyNoise(ceq_, noise);
                case FeatureName.RANDOM_NAN.value
                    % Similar to the case in the modifier_fun method.
                    rand_stream_random_nan = obj.payloadStream(seed, ceqCell{:}, xCell{:}, n_eval_ceq);
                    ceq_(rand_stream_random_nan.rand(size(ceq_)) < obj.options.(FeatureOptionKey.NAN_RATE.value)) = NaN;
                case FeatureName.TRUNCATED.value
                    % Similar to the case in the modifier_fun method.
                    rand_stream_truncated = obj.payloadStream(seed, ceqCell{:}, xCell{:}, n_eval_ceq);
                    digits = zeros(size(ceq_));
                    digits(ceq_ == 0) = obj.options.(FeatureOptionKey.SIGNIFICANT_DIGITS.value) - 1;
                    digits(ceq_ ~= 0) = obj.options.(FeatureOptionKey.SIGNIFICANT_DIGITS.value) - floor(log10(abs(ceq_(ceq_ ~= 0)))) - 1;
                    % Note that `floor(log10(abs(NaN))) = NaN` and `floor(log10(abs(Inf))) = Inf` in MATLAB.
                    for i_ceq = 1:numel(ceq_)
                        if ~isnan(ceq_(i_ceq)) && ~isinf(ceq_(i_ceq))
                            ceq_(i_ceq) = round(ceq_(i_ceq), digits(i_ceq));
                        end
                    end
                    if obj.options.(FeatureOptionKey.PERTURBED_TRAILING_DIGITS.value)
                        ceq_(ceq_ >= 0) = ceq_(ceq_ >= 0) + rand_stream_truncated.rand(size(ceq_(ceq_ >= 0))) .* (10 .^ (-digits(ceq_ >= 0)));
                        ceq_(ceq_ < 0) = ceq_(ceq_ < 0) - rand_stream_truncated.rand(size(ceq_(ceq_ < 0))) .* (10 .^ (-digits(ceq_ < 0)));
                    end
                case FeatureName.NONQUANTIFIABLE_CONSTRAINTS.value
                    % Set the elements whose absolute value are less than or equal to 10^(-6) to 0.
                    ceq_(abs(ceq_) <= 1e-6) = 0;
                    % Preserve NaN, matching Python and the raw truth channel.
                    ceq_(abs(ceq_) > 1e-6) = 1;
                case FeatureName.QUANTIZED.value
                    % Similar to the case in the modifier_fun method.
                    mesh_size = obj.options.(FeatureOptionKey.MESH_SIZE.value);
                    if strcmp(obj.options.(FeatureOptionKey.MESH_TYPE.value), 'relative')
                        mesh_size = mesh_size .* max(1, abs(x));
                    end
                    x = mesh_size .* round(x ./ mesh_size);
                    ceq_ = problem.ceq(x);
                otherwise
                    % Do nothing
            end
        end

        function noise = computeNoise(obj, x, seed, n_eval, base_values, noise_size)
            if strcmp(obj.options.(FeatureOptionKey.NOISE_MODE.value), 'deterministic')
                noise = obj.evaluateNoiseMap(x);
                noise = noise .* ones(noise_size);
                return;
            end

            xCell = num2cell(x);
            if iscell(base_values)
                baseCell = base_values;
            else
                baseCell = num2cell(base_values);
            end
            rand_stream_noisy = obj.payloadStream(seed, baseCell{:}, xCell{:}, n_eval);
            if strcmp(obj.options.(FeatureOptionKey.DISTRIBUTION.value), 'gaussian')
                noise = randn(rand_stream_noisy, noise_size);
            elseif strcmp(obj.options.(FeatureOptionKey.DISTRIBUTION.value), 'uniform')
                noise = 2 * rand(rand_stream_noisy, noise_size) - 1;
            else
                noise = obj.options.(FeatureOptionKey.DISTRIBUTION.value)(rand_stream_noisy, noise_size);
            end
        end

        function value = applyNoise(obj, value, noise)
            if strcmp(obj.options.(FeatureOptionKey.NOISE_TYPE.value), 'absolute')
                value = value + obj.options.(FeatureOptionKey.NOISE_LEVEL.value) * noise;
            elseif strcmp(obj.options.(FeatureOptionKey.NOISE_TYPE.value), 'relative')
                value = value .* (1.0 + obj.options.(FeatureOptionKey.NOISE_LEVEL.value) * noise);
            else
                % We need the noise to be symmetric with respect to 0.
                value = value + max(1, abs(value)) .* (obj.options.(FeatureOptionKey.NOISE_LEVEL.value) * noise);
            end
        end

        function noise = evaluateNoiseMap(obj, x)
            noise_map = obj.options.(FeatureOptionKey.NOISE_MAP.value);
            if isa(noise_map, 'function_handle')
                noise = noise_map(x);
            elseif strcmp(noise_map, 'chebyshev')
                noise = obj.chebyshevNoiseMap(x);
            else
                error("MATLAB:Feature:noise_map_InvalidInput", "Option `noise_map` must be 'chebyshev' or a function handle.")
            end
            if ~isrealscalar(noise)
                error("MATLAB:Feature:noise_map_InvalidOutput", "The output of `noise_map` must be a real scalar.")
            end
            noise = double(noise);
        end


        function rand_stream = payloadStream(obj, seed, varargin)
            % Per-query stream seeded by the run/stage seed and the observed
            % payload (values, point, served index). The legacy product mixer
            % is kept unchanged for the identity/single strategies; composed
            % views use the word fold of matlab-stage-horner32-v2.
            if strcmp(obj.payload_mixer, 'horner32-words')
                rand_stream = optiprofiler_internal.FeatureKernel.horner32_payload_rng(seed, varargin{:});
            else
                rand_stream = obj.default_rng(seed, varargin{:});
            end
        end
    end

    methods (Static)
        function rand_stream = horner32_payload_rng(seed, varargin)
            % Composed-views payload mixer (seed_policy matlab-stage-horner32-v2).
            % Starting from the 32-bit stage/channel seed, fold the IEEE-754
            % words (two uint32 per double, native byte order) of every payload
            % element in argument order with the exact 32-bit Horner rule
            % state = mod(65599 * state + word, 2^32), the same fold that
            % deriveFeatureStageSeed applies to stage identities. Each element
            % contributes at its own position, so a zero coordinate, a zero
            % value or a zero counter no longer removes the dependence on the
            % other elements as the legacy product mixer does. Negative zero
            % is folded as zero and every NaN as one canonical pattern. This is
            % a finite 32-bit hash: distinct payloads can collide, and it makes
            % no claim of statistical independence between streams.
            if ~(isnumeric(seed) && isreal(seed) && isscalar(seed)) || ~isfinite(seed)
                seed = 0;
            end
            state = mod(floor(double(seed)), 2^32);
            for k = 1:numel(varargin)
                values = double(varargin{k});
                values = reshape(values, 1, []);
                values(values == 0) = 0;
                values(isnan(values)) = NaN;
                words = double(typecast(values, 'uint32'));
                for word = words
                    state = mod(65599 * state + word, 2^32);
                end
            end
            rand_stream = RandStream('mt19937ar', 'Seed', state);
        end

        function noise = chebyshevNoiseMap(x)
            % Deterministic noise map from Moré and Wild, "Benchmarking
            % derivative-free optimization algorithms" (2009).
            alpha = 0.9 * sin(100 * norm(x, 1)) * cos(100 * norm(x, Inf)) + 0.1 * cos(norm(x, 2));
            noise = alpha .* (4 * alpha.^2 - 3);
        end

        function rand_stream = default_rng(seed, varargin)
            % Generate a random number generator.
            %
            % Parameters
            % ----------
            % seed : double, but default one is 'shuffle'
            %     Seed used to generate an initial random number generator.
            % varargin : array of double
            %     Arguments used to generate the returned random number generator.
            %
            % Returns
            % -------
            % rand_stream : RandStream
            %     Random number generator.

            % Create an initial rand_stream with the given seed
            if nargin < 1 || isempty(seed)
                seed = 'shuffle';  % Default behavior is to shuffle the seed.
            end
            if ~strcmp(seed, 'shuffle')
                if isnan(seed) || isinf(seed)
                    seed = 0;
                end
                if ~isrealscalar(seed)
                    error("MATLAB:Feature:SeedNotEvenReal", "The input `seed` for method `default_rng` in `Feature` should be a real number.")
                end
                seed = mod(floor(seed), 2^32);
            end

            % Convert all elements in varargin to double
            varargin = cellfun(@double, varargin, 'UniformOutput', false);

            % Create a RandStream object with the specified seed
            rand_stream = RandStream('mt19937ar', 'Seed', seed);

            % Generate a new seed based on the initial rand_stream and the additional arguments
            newSeed = abs(sin(1e5 * randn(rand_stream, 1)) + sum(sin(1e5 * prod(cellfun(@(x) x, varargin))))) * 1e9;
            newSeed = mod(floor(newSeed), 2^32);
            if isnan(newSeed) || isinf(newSeed)
                newSeed = 0;
            end
            rand_stream = RandStream('mt19937ar', 'Seed', floor(newSeed));
        end
    end
end

function tf = isrealscalar(x)
% Same predicate as src/private/isrealscalar, without package path coupling.
    tf = isnumeric(x) && isreal(x) && isscalar(x);
end

% Affine changes of variables x = A * y + b (linearly_transformed and custom
% with mod_affine).
%
% The bounds of the transformed problem have two representations. If A is
% diagonal they stay bounds, scaled by diag(inv) (the diagonal shortcut). For
% every invertible A they can be posed as linear rows of A (the generic
% representation), which needs A only. The bounds used to be classified by
% isdiag(inv) and the linear rows by isdiag(A). A diagonal A with one
% roundoff-sized off-diagonal entry in inv then made the bounds infinite while
% no bound row was added: the bounds left the posed problem without a word, and
% the truth went on scoring them. Hence:
%
% - one transformation per problem and seed, produced and validated once and
%   read by the initial point, the bounds, both kinds of linear constraints and
%   every evaluation (modifier_affine). User code may have a state, and a
%   second answer may be another map;
% - one decision, read by the bounds and by both kinds of linear constraints
%   (affineIsDiagonal). It is exact. The feasible set is a box exactly when A
%   is diagonal; an off-diagonal entry couples two variables, and what it moves
%   depends on the size of the other variable, which no tolerance on the entry
%   knows (4e-16 next to bounds of 1e16 moves the set by 4). inv cannot change
%   the set and decides nothing, except that the shortcut scales by diag(inv)
%   and so requires it to be the reciprocal of diag(A) to roundoff. Whatever
%   is not the shortcut is the generic representation, so the decision chooses
%   a representation and never whether a bound is posed;
% - everything the decision rests on is validated first, and what cannot be
%   represented raises instead of being approximated (checkedAffine, shifted,
%   scaledBounds, composedRows, pulledBack, refuseToReplaceBoundRows). Failing
%   closed matters here because the loss is silent: the solver is handed an
%   easier problem and is then scored on the original. Representable means
%   both ends of the range. Overflow: a finite quantity has to stay finite.
%   Underflow: a nonzero bound has to stay a normal number (at least realmin,
%   2.2e-308), and an entry of a transported row or right-hand side must not
%   have all of its terms below realmin; below it the spacing of numbers is
%   absolute, so digits are lost, and at zero an interval collapses to a point
%   and a row to no constraint. An entry that is zero because normal terms
%   cancel has lost nothing and is not refused;
% - the initial point is verified where it is used: no tolerance on the
%   matrices bounds an error at a point (pulledBack);
% - derivatives follow the same map. FeaturedProblem.grad and the other five
%   are those of the original callbacks in the variables of the solver, by the
%   chain rule of x = A * y + b; without a change of variables they are the
%   established passthrough, and a composition provides none.

function [A, b, inv] = checkedAffine(A, b, inv, n, supplied)
% The change of variables x = A * y + b as validated full double arrays.
% Checked in this order, each with its own message: real arrays of sizes
% n-by-n, n and n-by-n; finite entries (NaN fails no inequality, so it has to
% be asked for); numerical invertibility,
% norm(abs(inv) * abs(A), inf) < 1 / eps; and, if inv was SUPPLIED by user
% code, consistency from both sides (see below).
%
% This condition number does not change when the rows of A (the units of the
% original variables) are scaled, and 1 / eps is where a matrix is singular to
% working precision. The usual norm(A) * norm(inv) would refuse exact
% transformations that are merely badly scaled: here a diagonal scaling has
% condition number 1 whatever its entries, and the scaled rotation of
% linearly_transformed at most n. The framework's own inverse of a rotation is
% exact up to roundoff times the usual condition number and is not held to the
% residual test.
%
% Consistency: norm(E, 'fro') <= 1e-8 * n for
% E = max(abs(A * inv - I) - 64 * n * eps * abs(A) * abs(inv), 0) and for
% E = max(abs(inv * A - I) - 64 * n * eps * abs(inv) * abs(A), 0). Both products are needed:
% with A = diag(1e-8, 1e8) the first is an identity to 1e-16 for an inv with
% which the second misses it by 1. Each entry is measured against the terms it
% is summed from, which is what a change of units scales: a rotation with one
% variable in units of 1e13 misses the identity by 1e-3 in one of the products,
% new variable or original one, and is consistent to roundoff. The allowance is
% the rounding level of the products themselves (see roundingAllowance), not a
% relative 1e-8: an inverse of an ill-conditioned matrix with one entry off by
% 1e-9 of its size is refused, as it was by the plain
% norm(A * inv - I, 'fro') <= 1e-8 * n, which this rule is never stricter than.
% No rule on the matrices bounds an error at a point, so the one point that inv
% is used for is verified there (pulledBack).
%
% Identifiers: AffineTransformationInvalid when A or b is not usable data;
% AffineTransformationNotInvertible whenever inv cannot be the inverse of A
% (wrong size or type, not finite, numerically singular, inconsistent).
    usable = @(value) (isnumeric(value) || islogical(value)) && isreal(value);
    if ~(usable(A) && isequal(size(A), [n, n]))
        error("MATLAB:Feature:AffineTransformationInvalid", "The affine transformation matrix must be a real matrix of size %d-by-%d.", n, n);
    end
    if ~(usable(b) && numel(b) == n && (isvector(b) || n == 0))
        error("MATLAB:Feature:AffineTransformationInvalid", "The affine transformation vector must be a real vector of size %d.", n);
    end
    if ~(usable(inv) && isequal(size(inv), [n, n]))
        error("MATLAB:Feature:AffineTransformationNotInvertible", "The inverse of the affine transformation matrix must be a real matrix of size %d-by-%d.", n, n);
    end
    A = full(double(A));
    b = full(double(b(:)));
    inv = full(double(inv));
    if ~(all(isfinite(A(:))) && all(isfinite(b)))
        error("MATLAB:Feature:AffineTransformationInvalid", "The affine transformation matrix and vector must be finite.");
    end
    if ~all(isfinite(inv(:)))
        error("MATLAB:Feature:AffineTransformationNotInvertible", "The inverse of the affine transformation matrix must be finite.");
    end
    % Finite entries can still overflow in a product. All tests are written so
    % that an infinite or NaN result fails them.
    terms = abs(inv) * abs(A);
    condition = norm(terms, inf);
    if ~(condition * eps < 1)
        error("MATLAB:Feature:AffineTransformationNotInvertible", "The affine transformation is numerically singular: norm(abs(inv) * abs(A), inf) is %.3g, not below 1 / eps.", condition);
    end
    if supplied && ~(identityResidual(A * inv, abs(A) * abs(inv), n) <= 1e-8 * n && identityResidual(inv * A, terms, n) <= 1e-8 * n)
        error("MATLAB:Feature:AffineTransformationNotInvertible", "The multiplication of the affine transformation matrix and its inverse is not an identity matrix.");
    end
end

function value = identityResidual(product, terms, n)
% What rounding cannot explain: each entry is forgiven the rounding level of
% the products it is summed from, and nothing more. NaN if the terms overflowed.
    if ~all(isfinite(terms(:)))
        value = NaN;
        return;
    end
    value = norm(max(abs(product - eye(n)) - roundingAllowance() * n * eps * terms, 0), 'fro');
end

function tf = affineIsDiagonal(A, inv)
% The one structural decision: whether the bounds stay bounds.
% They do if A is diagonal, exactly, and diag(inv) is the reciprocal of diag(A)
% to roundoff, abs(inv(i, i) * A(i, i) - 1) <= 8 * eps. The first is what makes
% the feasible set a box. The second is what makes diag(inv) times a bound that
% bound in the new variables: two factors that are each rounded to within two
% units in the last place differ from exact reciprocals by less than 5 * eps
% (the pair that linearly_transformed builds is within 1.5 * eps), whereas a
% supplied inverse is only held to 1e-8, and a bound may be of any size.
% Off-diagonal entries of inv change no feasible set and are not read.
% Everything else takes the generic representation, which poses the bounds
% exactly whatever inv is.
    tf = isdiag(A) && all(abs(diag(inv) .* diag(A) - 1) <= 8 * eps);
end

function refuseToReplaceBoundRows(kernel, seed, problem, key, fixed)
% A supplied linear modifier replaces the linear constraints verbatim. Under a
% custom affine map that is not diagonal, and without mod_bounds, the framework
% poses the bounds as exactly those constraints (FIXED variables as
% equalities, the other finite bounds as inequalities), so replacing them would
% drop the bounds silently. There is then no representation left for them.
    if ~isfield(kernel.options, FeatureOptionKey.MOD_AFFINE.value) || isfield(kernel.options, FeatureOptionKey.MOD_BOUNDS.value)
        return;
    end
    [A, ~, inv] = kernel.modifier_affine(seed, problem);
    if affineIsDiagonal(A, inv)
        return;  % the bounds stay bounds
    end
    is_fixed = problem.xl == problem.xu;
    if fixed
        at_stake = is_fixed;
    else
        at_stake = (isfinite(problem.xl) | isfinite(problem.xu)) & ~is_fixed;
    end
    if any(at_stake)
        error("MATLAB:Feature:AffineBoundsNotRepresentable", "The affine transformation is not diagonal (or the diagonal of its inverse is not the reciprocal of its diagonal to roundoff), so the finite bounds of the problem are posed as linear constraints, which %s replaces: the bounds would be dropped silently. Supply mod_bounds as well, or a diagonal transformation with its exact inverse.", key);
    end
end

function values = shifted(values, b)
% values - b for bounds, where an infinite entry means "none"; a finite entry
% has to stay finite.
    was_finite = isfinite(values);
    values = values - b;
    if any(was_finite & ~isfinite(values))
        error("MATLAB:Feature:AffineBoundsNotRepresentable", "A finite bound is not representable after the affine transformation: it overflows, and posing it as infinite would drop it silently.");
    end
end

function [xl, xu] = scaledBounds(scale, lower, upper)
% Bounds of the diagonal shortcut, scale .* [lower, upper], swapped where the
% scale is negative. A finite bound has to stay finite, and a nonzero one a
% normal number: 1e-200 * [1e-200, 2e-200] is [0, 0].
    scaled_lower = scale .* lower;
    scaled_upper = scale .* upper;
    % A finite bound times a finite scale can overflow. It would then be posed
    % as "no bound", which is the silent loss this file rules out.
    if any(isfinite(lower) & ~isfinite(scaled_lower)) || any(isfinite(upper) & ~isfinite(scaled_upper))
        error("MATLAB:Feature:AffineBoundsNotRepresentable", "A finite bound is not representable after the affine transformation: it overflows, and posing it as infinite would drop it silently.");
    end
    % And it can underflow: a nonzero bound below the smallest normal number
    % has lost its digits or is zero, and an interval can collapse to a point.
    if any(lower ~= 0 & abs(scaled_lower) < realmin) || any(upper ~= 0 & abs(scaled_upper) < realmin)
        error("MATLAB:Feature:AffineBoundsNotRepresentable", "A finite bound is not representable after the affine transformation: it underflows (its magnitude falls below the smallest normal number), and posing it as zero would move it.");
    end
    xl = min(scaled_lower, scaled_upper);
    xu = max(scaled_lower, scaled_upper);
end

function [rows, moved] = composedRows(matrix, rhs, A, b)
% The linear constraints with coefficients MATRIX and right-hand side RHS in the
% new variables: matrix * A and rhs - matrix * b (rhs itself if b is empty:
% there is no shift). A row of finite data has to stay finite: an infinite
% right-hand side is counted as no constraint, a negative one is satisfied
% nowhere, and NaN compares as satisfied. And no entry may be lost to
% underflow: a coefficient or a right-hand side that has a product of nonzero
% factors among its terms while the sum of the absolute terms is below the
% smallest normal number (the right-hand side itself counts as a term of the
% shifted one).
    rows = matrix * A;
    moved = rhs;
    if ~isempty(b)
        moved = rhs - matrix * b;
    end
    was_finite = all(isfinite(matrix), 2) & isfinite(rhs);
    if ~(all(all(isfinite(rows(was_finite, :)))) && all(isfinite(moved(was_finite))))
        error("MATLAB:Feature:AffineLinearConstraintsNotRepresentable", "A linear constraint is not representable after the affine transformation: it overflows, and a coefficient or a right-hand side that is not finite is not the constraint.");
    end
    % Underflow: an entry that has a product of nonzero factors among its terms
    % while the sum of the absolute terms is below the smallest normal number
    % is lost (a row of zeros is no constraint). Cancellation is not underflow:
    % it leaves the absolute terms normal.
    M = matrix(was_finite, :);
    lost = (double(M ~= 0) * double(A ~= 0) > 0) & (abs(M) * abs(A) < realmin);
    lost_rhs = false;
    if ~isempty(b)
        lost_rhs = (double(M ~= 0) * double(b ~= 0) > 0) & (abs(rhs(was_finite)) + abs(M) * abs(b) < realmin);
    end
    if any(lost(:)) || any(lost_rhs(:))
        error("MATLAB:Feature:AffineLinearConstraintsNotRepresentable", "A linear constraint is not representable after the affine transformation: a coefficient or a right-hand side underflows (all of its terms fall below the smallest normal number).");
    end
end

function point = pulledBack(A, inv, x0, b)
% The initial point in the new variables, and the proof that it is one.
% inv * (x0 - b) is the point if it is mapped back to x0 to the rounding of that
% evaluation, in every component:
% abs(A * y + b - x0) <= 64 * n * eps * (abs(A) * abs(y) + abs(b) + abs(x0)).
% No tolerance on the matrices can stand in for this test: A = I with
% inv(1, 2) = 1e-12 is an identity to 1e-12 from both sides and moves
% x0 = (0, 1e14) to (100, 1e14), 99 outside bounds of [-1, 1]. If the point
% fails the test, the equation A * y = x0 - b is solved instead, which needs A
% only; if that point fails it as well (the shift or a product overflows, a
% component underflows, or A is too ill conditioned at this point),
% construction raises. A point that passes is kept bitwise, so an inverse that
% is good to roundoff at x0 (the framework's own, measured at most 3.1 units of
% the allowance; or a supplied one) gives the point it always gave. The test is
% by component on purpose: a norm would let a component of 1e14 excuse an error
% of 100 in another one. b is empty if there is no shift.
    n = numel(x0);
    shift = zeros(n, 1);
    point = inv * x0;
    if ~isempty(b)
        shift = b;
        point = inv * (x0 - b);
    end
    mapped = @(y) all(abs(A * y + shift - x0) <= roundingAllowance() * n * eps * (abs(A) * abs(y) + abs(shift) + abs(x0)));  % false for NaN and for an infinite component
    if ~all(isfinite(x0)) || mapped(point)
        return;  % (a point that is not finite is the user's: it is transported as it is)
    end
    states = [warning('off', 'MATLAB:nearlySingularMatrix'), warning('off', 'MATLAB:singularMatrix')];
    restore = onCleanup(@() warning(states));
    point = A \ (x0 - shift);
    if ~all(isfinite(point))
        error("MATLAB:Feature:AffineInitialPointNotRepresentable", "The initial point is not representable after the affine transformation: it overflows.");
    end
    if ~mapped(point)
        error("MATLAB:Feature:AffineInitialPointNotRepresentable", "The initial point is not representable after the affine transformation: the point found in the new variables is not mapped back to it to roundoff (it underflows, or the transformation is too ill conditioned at that point).");
    end
end

function value = roundingAllowance()
% What counts as the rounding of a sum of products: this many times n * eps
% times the sum of the absolute values of its terms. Measured over pairs that
% are consistent to roundoff (built by formula as linearly_transformed builds
% its own, or by LU, condition numbers up to 1e14, NumPy 1.24 to 2.5 and MATLAB
% R2026a): at most 6.6 for the products of checkedAffine and 3.1 for the point
% of pulledBack. An inverse that is off by 1e-9 of an entry is at 4.5e6.
    value = 64;
end
