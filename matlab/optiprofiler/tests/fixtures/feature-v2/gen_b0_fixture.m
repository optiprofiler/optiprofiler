function gen_b0_fixture(output_root)
%GEN_B0_FIXTURE writes the frozen 1.x FeaturedProblem fixture of this folder.
%
%   It has to be run with the 1.x implementation on the path and nothing of a
%   later one (see b0-featured-problem.README.md for the exact commands): the
%   fixture is evidence of what that implementation saved, so it is never to be
%   regenerated with newer code.
%
%   The two callbacks are built with str2func on purpose. An anonymous function
%   written in a file records the absolute path of that file, and loading it on
%   a machine without that path warns (MATLAB:dispatcher:UnresolvedFunctionHandle).
%   The first capture of this fixture recorded the path of this generator on the
%   machine that ran it, so the tests that require a silent load passed there
%   and failed everywhere else. A handle from str2func records no file.

    if ~isfolder(output_root), mkdir(output_root); end
    fprintf('Feature from: %s\n', which('Feature'));
    fun = str2func('@(x) sum(x.^2)');
    cub = str2func('@(x) x(1) - 0.5');
    for handle = {fun, cub}
        info = functions(handle{1});
        assert(isempty(info.file), 'gen_b0_fixture:CreatorRecorded', 'A callback records the file %s.', info.file);
    end
    problem = Problem(struct('fun', fun, 'x0', [1; 2], 'cub', cub, 'name', 'b0fixture'));
    feature_noisy = Feature('noisy', 'noise_level', 0.5, 'n_runs', 2);
    fp_noisy = FeaturedProblem(problem, feature_noisy, 6, 17);
    observed = struct();
    observed.noisy_f1 = fp_noisy.fun([1; 2]);
    observed.noisy_f2 = fp_noisy.fun([0.5; 1]);
    observed.noisy_c1 = fp_noisy.cub([1; 2]);
    fp_plain = FeaturedProblem(problem, Feature('plain'), 6, 17);
    observed.plain_f1 = fp_plain.fun([1; 2]);
    save(fullfile(output_root, 'b0-featured-problem.mat'), 'fp_noisy', 'fp_plain', '-v7');
    % Continue the same live objects under 1.x: the values a restored object must reproduce.
    continuation = struct();
    continuation.noisy_f3 = fp_noisy.fun([0.25; 0.5]);
    continuation.noisy_c2 = fp_noisy.cub([0.25; 0.5]);
    continuation.noisy_maxcv = fp_noisy.maxcv([0.25; 0.5]);
    continuation.noisy_fun_hist = fp_noisy.fun_hist;
    continuation.noisy_cub_hist = fp_noisy.cub_hist;
    continuation.noisy_n_eval_fun = fp_noisy.n_eval_fun;
    continuation.plain_f2 = fp_plain.fun([0.25; 0.5]);
    continuation.plain_fun_hist = fp_plain.fun_hist;
    continuation.observed = observed;
    continuation.matlab = version;
    continuation.source_sha = '1c4b7d9170ad40e1112e33e5a5638cecc9e644eb';
    save(fullfile(output_root, 'b0-featured-problem-continuation.mat'), 'continuation', '-v7');
    disp(continuation);
    fprintf('B0-FIXTURE-WRITTEN\n');
end
