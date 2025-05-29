function success = setupPool(n_jobs)
%SETUPPOOL set up the parallel pool by the number of jobs.

    % If `n_jobs` is 1, then do nothing.
    if n_jobs == 1
        fprintf('\nINFO: The option `n_jobs` is set to 1. We will not use parallel computing.\n');
        success = false;
        return;
    end

    % Check whether the Parallel Computing Toolbox is available.
    if exist('parpool', 'file') ~= 2
        fprintf('\nINFO: Parallel Computing Toolbox is not available. We will not use parallel computing.\n');
        success = false;
        return;
    end

    % Check whether there is an existing parallel pool.
    if ~isempty(gcp('nocreate'))
        myCluster = gcp('nocreate');
        n_workers = myCluster.NumWorkers;
        if n_jobs <= n_workers
            % If there are more workers than `n_jobs`, then we will use `n_jobs` workers.
            fprintf('\nINFO: There is an existing parallel pool with %d workers. We will use %d workers by the option `n_jobs`.\n', n_workers, n_jobs);
            success = true;
            return;
        else
            % If there are fewer workers than `n_jobs`, then we will close the existing pool and try to open a new one with maximum workers later.
            fprintf('\nINFO: There is an existing parallel pool with %d workers, which is fewer than the option `n_jobs` (%d). We will close the existing pool and try to open a new one with maximum workers.\n', n_workers, n_jobs);
            delete(myCluster);
        end
    end

    % Try to open a new parallel pool with maximum workers.
    try
        myCluster = parcluster();
        n_workers = myCluster.NumWorkers;
        if n_jobs > n_workers
            fprintf('\nINFO: The option `n_jobs` (%d) is greater than the maximum number of workers (%d) in your cluster profile setting. We will use %d workers instead.\n', n_jobs, n_workers, n_workers);
            fprintf("\nINFO: You may change the maximum number of workers in your cluster by running following command in the MATLAB command window:\n\n");
            fprintf("    myCluster = parcluster(<name_of_your_cluster>); myCluster.NumWorkers = <new_number_of_workers>; saveProfile(myCluster);\n");
        end
        n_jobs = min(n_jobs, n_workers);
        parpool(n_jobs);
        success = true;
    catch ME
        fprintf('\nINFO: Failed to open a parallel pool with %d workers. Error message: %s\n', n_jobs, ME.message);
        success = false;
    end
end