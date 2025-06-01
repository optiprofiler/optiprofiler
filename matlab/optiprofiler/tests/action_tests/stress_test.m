function stress_test(benchmark_id)
%STRESS_TEST tests randomly high-dimensional (unconstrained) problems from s2mpj or matcutest.

    % Go to the directory of this repository.
    cd(fullfile(fileparts(mfilename('fullpath')), '../../../../output'));

    % Get the list of problems from s2mpj and matcutest whose dimension ranging from 2000 to 10000.
    s2mpj_pb_list = s2mpj_select(struct('ptype', 'u', 'mindim', 2000, 'maxdim', 10000));
    if isunix && ~ismac
        matcutest_pb_list = matcutest_select(struct('ptype', 'u', 'mindim', 2000, 'maxdim', 10000));
        pb_list = [s2mpj_pb_list matcutest_pb_list];
    else
        pb_list = s2mpj_pb_list;
    end

    % Use the current wall-clock time to randomly select 5 problems.
    time = datetime;
    seed = 100*mod(year(time), 100) + week(time);
    fprintf("Seed: %d\n\n", seed);
    rand_stream = RandStream("mt19937ar", "Seed", seed);
    random_idx = rand_stream.randperm(numel(pb_list), 5);
    selected_pb_list = pb_list(random_idx);
    fprintf("Selected problem names: %s\n\n", strjoin(selected_pb_list, ', '));

    % Set options for OptiProfiler.
    options.seed = seed;
    options.problem_names = selected_pb_list;
    options.mindim = 2000;
    options.maxdim = 10000;
    options.max_eval_factor = 0.01;
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