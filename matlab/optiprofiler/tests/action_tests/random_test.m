function random_test(benchmark_id)
%RANDOM_TEST tests OptiProfiler with random options.

    % Use the current wall-clock time to generate a random seed.
    time = datetime;
    seed = 100*mod(year(time), 100) + month(time) + week(time) + day(time);
    fprintf("Seed: %d\n\n", seed);
    rand_stream = RandStream("mt19937ar", "Seed", seed);

    % Use the random seed to generate random options.
    options = struct();

    n_jobs_choices = (1:5);
    options.n_jobs = n_jobs_choices(rand_stream.randi(length(n_jobs_choices), 1, 1));

    options.seed = seed;

    options.benchmark_id = benchmark_id;

    solvers_choices = {@fmincon_test1, @fmincon_test2, @fmincon_test3};
    solver_names_choices = {'sqp', 'interior-point', 'active-set'};
    solver_2_or_3 = rand_stream.randi([2, 3], 1, 1);
    if solver_2_or_3 == 2
        idx = rand_stream.randperm(2, 2);
        solvers = solvers_choices(idx);
        options.solver_names = solver_names_choices(idx);
    else
        solvers = solvers_choices;
        options.solver_names = solver_names_choices;
    end

    error_bar_choices = {'minmax', 'meanstd'};
    options.errorbar_type = error_bar_choices{rand_stream.randi(length(error_bar_choices), 1, 1)};

    max_tol_order_choices = (1:16);
    options.max_tol_order = max_tol_order_choices(rand_stream.randi(length(max_tol_order_choices), 1, 1));

    options.max_eval_factor = rand_stream.rand(1, 1) * 1000;

    options.project_x0 = (rand_stream.rand(1, 1) < 0.5);

    options.run_plain = (rand_stream.rand(1, 1) < 0.5);

    options.score_only = (rand_stream.rand(1, 1) < 0.5);

    options.summarize_performance_profiles = (rand_stream.rand(1, 1) < 0.5);

    options.summarize_data_profiles = (rand_stream.rand(1, 1) < 0.5);

    if length(solvers) == 2
        options.summarize_log_ratio_profiles = (rand_stream.rand(1, 1) < 0.5);
    else
        options.summarize_log_ratio_profiles = false;
    end

    options.summarize_output_based_profiles = (rand_stream.rand(1, 1) < 0.5);

    options.solver_verbose = rand_stream.randi([0, 2], 1, 1);

    options.semilogx = (rand_stream.rand(1, 1) < 0.5);

    options.normalized_scores = (rand_stream.rand(1, 1) < 0.5);

    line_colors_choices = {'r', 'g', 'b', 'c', 'm', 'y', 'k'};
    options.line_colors = line_colors_choices(rand_stream.randperm(length(line_colors_choices), length(solvers)));

    line_styles_choices = {...
    '-', '-o', '-+', '-*', '-.', '-x', '-s', '-d', '-^', '-v', '->', '-<', '-p', '-h', 'o-', '+-', '*-', '.-', 'x-', 's-', 'd-', '^-', 'v-', '>-', '<-', 'p-', 'h-',...
    '-.', '-.o', '-.+', '-.*', '-..', '-.x', '-.s', '-.d', '-.^', '-.v', '-.>', '-.<', '-.p', '-.h', 'o-.', '+-.', '*-.', '.-.', 'x-.', 's-.', 'd-.', '^-.', 'v-.', '>-.', '<-.', 'p-.', 'h-.',...
    ':', ':o', ':+', ':*', ':.', ':x', ':s', ':d', ':^', ':v', ':>', ':<', ':p', ':h', 'o:', '+:', '*:', '.:', 'x:', 's:', 'd:', '^:', 'v:', '>:', '<:', 'p:', 'h:',...
    '--', '--o', '--+', '--*', '--.', '--x', '--s', '--d', '--^', '--v', '-->', '--<', '--p', '--h', 'o--', '+--', '*--', '.--', 'x--', 's--', 'd--', '^--', 'v--', '>--', '<--', 'p--', 'h--'...
    };
    options.line_styles = line_styles_choices(rand_stream.randperm(length(line_styles_choices), length(solvers)));

    options.line_widths = rand_stream.rand(1, length(solvers)) * 2 + 0.5;

    bar_colors_choices = {'r', 'g', 'b', 'c', 'm', 'y', 'k'};
    options.bar_colors = bar_colors_choices(rand_stream.randperm(length(bar_colors_choices), length(solvers)));

    feature_name_choices = {'plain', 'perturbed_x0', 'noisy', 'truncated', 'permuted', 'linearly_transformed', 'random_nan', 'unrelaxable_constraints', 'nonquantifiable_constraints', 'quantized'};
    options.feature_name = feature_name_choices{rand_stream.randi(length(feature_name_choices), 1, 1)};

    switch options.feature_name
        case 'plain'
            options.n_runs = 1;
        case 'perturbed_x0'
            options.n_runs = rand_stream.randi([1, 3], 1, 1);
            options.perturbation_level = rand_stream.rand(1, 1);
            distribution_choices = {'gaussian', 'spherical'};
            options.distribution = distribution_choices{rand_stream.randi(length(distribution_choices), 1, 1)};
        case 'noisy'
            options.n_runs = rand_stream.randi([1, 3], 1, 1);
            options.noise_level = rand_stream.rand(1, 1);
            distribution_choices = {'gaussian', 'uniform'};
            options.distribution = distribution_choices{rand_stream.randi(length(distribution_choices), 1, 1)};
            noise_type_choices = {'absolute', 'relative', 'mixed'};
            options.noise_type = noise_type_choices{rand_stream.randi(length(noise_type_choices), 1, 1)};
        case 'truncated'
            options.n_runs = rand_stream.randi([1, 3], 1, 1);
            options.significant_digits = rand_stream.randi([1, 10], 1, 1);
            options.perturbed_trailing_digits = (rand_stream.rand(1, 1) < 0.5);
        case 'linearly_transformed'
            options.n_runs = rand_stream.randi([1, 3], 1, 1);
            options.rotated = (rand_stream.rand(1, 1) < 0.5);
            options.condition_factor = rand_stream.rand(1, 1) * 100;
        case 'random_nan'
            options.n_runs = rand_stream.randi([1, 3], 1, 1);
            options.nan_rate = rand_stream.rand(1, 1) * 0.5;
        case 'unrelaxable_constraints'
            options.n_runs = 1;
            options.unrelaxable_bounds = (rand_stream.rand(1, 1) < 0.5);
            options.unrelaxable_linear_constraints = (rand_stream.rand(1, 1) < 0.5);
            options.unrelaxable_nonlinear_constraints = (rand_stream.rand(1, 1) < 0.5);
        case 'nonquantifiable_constraints'
            options.n_runs = 1;
        case 'quantized'
            options.n_runs = 1;
            options.mesh_size = rand_stream.rand(1, 1);
            mesh_type_choices = {'absolute', 'relative'};
            options.mesh_type = mesh_type_choices{rand_stream.randi(length(mesh_type_choices), 1, 1)};
            options.ground_truth = (rand_stream.rand(1, 1) < 0.5);
    end

    if isunix() && ~ismac()
        options.plibs = {'s2mpj', 'matcutest'};
    else
        options.plibs = {'s2mpj'};
    end

    ptype_choices = {'u', 'b', 'l', 'n', 'ub', 'ul', 'un', 'bl', 'bn', 'ln', 'ubl', 'ubn', 'uln', 'bln', 'ubln'};
    options.ptype = ptype_choices{rand_stream.randi(length(ptype_choices), 1, 1)};

    options.mindim = rand_stream.randi([1, 3], 1, 1);
    options.maxdim = options.mindim + rand_stream.randi([1, 3], 1, 1);

    options.minb = 0;
    options.maxb = rand_stream.randi([20, 40], 1, 1);

    options.minlcon = 0;
    options.maxlcon = rand_stream.randi([10, 20], 1, 1);

    options.minnlcon = 0;
    options.maxnlcon = rand_stream.randi([10, 20], 1, 1);

    options.mincon = 0;
    options.maxcon = rand_stream.randi([20, 40], 1, 1);

    if rand_stream.rand(1, 1) < 0.5
        options.xlabel_data_profile = '';
    end

    if rand_stream.rand(1, 1) < 0.5
        options.ylabel_data_profile = '';
    end

    if rand_stream.rand(1, 1) < 0.5
        options.xlabel_performance_profile = '';
    end

    if rand_stream.rand(1, 1) < 0.5
        options.ylabel_performance_profile = '';
    end

    if rand_stream.rand(1, 1) < 0.5
        options.xlabel_log_ratio_profile = '';
    end

    if rand_stream.rand(1, 1) < 0.5
        options.ylabel_log_ratio_profile = '';
    end

    % Run the stress test.
    benchmark(solvers, options)
end

function x = fmincon_test1(varargin)
    options = optimoptions('fmincon', 'SpecifyObjectiveGradient', false, 'Algorithm', 'sqp');
    if nargin == 2
        fun = varargin{1};
        x0 = varargin{2};
        x = fmincon(fun, x0, [], [], [], [], [], [], [], options);
    elseif nargin == 4
        fun = varargin{1};
        x0 = varargin{2};
        xl = varargin{3};
        xu = varargin{4};
        x = fmincon(fun, x0, [], [], [], [], xl, xu, [], options);
    elseif nargin == 8
        fun = varargin{1};
        x0 = varargin{2};
        xl = varargin{3};
        xu = varargin{4};
        aub = varargin{5};
        bub = varargin{6};
        aeq = varargin{7};
        beq = varargin{8};
        x = fmincon(fun, x0, aub, bub, aeq, beq, xl, xu, [], options);
    elseif nargin == 10
        fun = varargin{1};
        x0 = varargin{2};
        xl = varargin{3};
        xu = varargin{4};
        aub = varargin{5};
        bub = varargin{6};
        aeq = varargin{7};
        beq = varargin{8};
        cub = varargin{9};
        ceq = varargin{10};
        nonlcon = @(x) deal(cub(x), ceq(x));
        x = fmincon(fun, x0, aub, bub, aeq, beq, xl, xu, nonlcon, options);
    end
end

function x = fmincon_test2(varargin)

    options = optimoptions('fmincon', 'SpecifyObjectiveGradient', false, 'Algorithm', 'interior-point');
    if nargin == 2
        fun = varargin{1};
        x0 = varargin{2};
        x = fmincon(fun, x0, [], [], [], [], [], [], [], options);
    elseif nargin == 4
        fun = varargin{1};
        x0 = varargin{2};
        xl = varargin{3};
        xu = varargin{4};
        x = fmincon(fun, x0, [], [], [], [], xl, xu, [], options);
    elseif nargin == 8
        fun = varargin{1};
        x0 = varargin{2};
        xl = varargin{3};
        xu = varargin{4};
        aub = varargin{5};
        bub = varargin{6};
        aeq = varargin{7};
        beq = varargin{8};
        x = fmincon(fun, x0, aub, bub, aeq, beq, xl, xu, [], options);
    elseif nargin == 10
        fun = varargin{1};
        x0 = varargin{2};
        xl = varargin{3};
        xu = varargin{4};
        aub = varargin{5};
        bub = varargin{6};
        aeq = varargin{7};
        beq = varargin{8};
        cub = varargin{9};
        ceq = varargin{10};
        nonlcon = @(x) deal(cub(x), ceq(x));
        x = fmincon(fun, x0, aub, bub, aeq, beq, xl, xu, nonlcon, options);
    end
end

function x = fmincon_test3(varargin)
    options = optimoptions('fmincon', 'SpecifyObjectiveGradient', false, 'Algorithm', 'active-set');
    if nargin == 2
        fun = varargin{1};
        x0 = varargin{2};
        x = fmincon(fun, x0, [], [], [], [], [], [], [], options);
    elseif nargin == 4
        fun = varargin{1};
        x0 = varargin{2};
        xl = varargin{3};
        xu = varargin{4};
        x = fmincon(fun, x0, [], [], [], [], xl, xu, [], options);
    elseif nargin == 8
        fun = varargin{1};
        x0 = varargin{2};
        xl = varargin{3};
        xu = varargin{4};
        aub = varargin{5};
        bub = varargin{6};
        aeq = varargin{7};
        beq = varargin{8};
        x = fmincon(fun, x0, aub, bub, aeq, beq, xl, xu, [], options);
    elseif nargin == 10
        fun = varargin{1};
        x0 = varargin{2};
        xl = varargin{3};
        xu = varargin{4};
        aub = varargin{5};
        bub = varargin{6};
        aeq = varargin{7};
        beq = varargin{8};
        cub = varargin{9};
        ceq = varargin{10};
        nonlcon = @(x) deal(cub(x), ceq(x));
        x = fmincon(fun, x0, aub, bub, aeq, beq, xl, xu, nonlcon, options);
    end
end