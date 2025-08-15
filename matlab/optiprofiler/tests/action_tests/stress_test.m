function stress_test(benchmark_id)
%STRESS_TEST tests randomly high-dimensional (unconstrained) problems from s2mpj or matcutest.

    % Get the list of problems from s2mpj and matcutest whose dimension ranging from 2000 to 10000.
    s2mpj_pb_list = s2mpj_select(struct('ptype', 'u', 'mindim', 2000, 'maxdim', 10000));
    if isunix && ~ismac
        matcutest_pb_list = matcutest_select(struct('ptype', 'u', 'mindim', 2000, 'maxdim', 10000));
        pb_list = [s2mpj_pb_list matcutest_pb_list];
    else
        pb_list = s2mpj_pb_list;
    end

    % Use the current wall-clock time to randomly select 5 problems.
    time_zone = 'Asia/Shanghai';
    dt = datetime('now', 'TimeZone', time_zone);
    seed = 100*mod(year(dt), 100) + week(dt);
    fprintf("Seed: %d\n\n", seed);
    rand_stream = RandStream("mt19937ar", "Seed", seed);
    random_idx = rand_stream.randperm(numel(pb_list), 3);
    selected_pb_list = pb_list(random_idx);
    fprintf("Selected problem names: %s\n\n", strjoin(selected_pb_list, ', '));

    % Set options for OptiProfiler.
    options.seed = seed;
    options.problem_names = selected_pb_list;
    options.mindim = 2000;
    options.maxdim = 10000;
    if isunix && ~ismac
        options.max_eval_factor = 0.8;
    else
        options.max_eval_factor = 0.2;
    end
    options.benchmark_id = benchmark_id;
    if isunix && ~ismac
        options.plibs = {'s2mpj', 'matcutest'};
    else
        options.plibs = {'s2mpj'};
    end
    solvers = {@fminunc, @fminsearch};

    % Run the stress test.
    benchmark(solvers, options)
end