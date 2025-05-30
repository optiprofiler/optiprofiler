function success = setupPool(n_jobs, silent)
%SETUPPOOL set up the parallel pool by the number of jobs.

    % Check whether the Parallel Computing Toolbox is available.
    if exist('parpool', 'file') ~= 2
        if ~silent
            fprintf('\nINFO: Parallel Computing Toolbox is not available. We will not use parallel computing.\n');
        end
        success = false;
        return;
    end

    % Check the default setting for the maximum number of workers in the cluster profile.
    defaultCluster = parcluster();
    max_workers = defaultCluster.NumWorkers;

    % Check whether there is an existing parallel pool.
    if ~isempty(gcp('nocreate'))
        myCluster = gcp('nocreate');
        n_workers = myCluster.NumWorkers;
        if n_jobs <= n_workers
            % If there are more workers than `n_jobs`, then we will use `n_jobs` workers.
            if ~silent
                fprintf('\nINFO: There is an existing parallel pool with %d workers. We will use %d workers by the option `n_jobs`.\n', n_workers, n_jobs);
            end
            success = true;
            return;
        elseif n_workers == max_workers
            % If there are maximum workers in the cluster profile, then we will use all workers.
            if ~silent
                fprintf('\nINFO: There is an existing parallel pool with %d workers, which is equal to the maximum number of workers in your cluster profile setting. We will use all %d workers instead of the option `n_jobs` (%d).\n', n_workers, n_workers, n_jobs);
                fprintf("\nINFO: You may change the maximum number of workers in your cluster by running following command in the MATLAB command window:\n\n");
                fprintf("    myCluster = parcluster(); myCluster.NumWorkers = <new_number_of_workers>; saveProfile(myCluster);\n\n");
            end
            success = true;
            return;
        else
            % If there are fewer workers than `n_jobs`, then we will close the existing pool and try to open a new one with more workers later.
            if ~silent
                fprintf('\nINFO: There is an existing parallel pool with %d workers, which is fewer than the option `n_jobs` (%d). We will close the existing pool and try to open a new one with more workers.\n', n_workers, n_jobs);
            end
            delete(myCluster);
        end
    end

    % Try to open a new parallel pool with maximum workers.
    try
        if n_jobs > max_workers
            if ~silent
                fprintf('\nINFO: The option `n_jobs` (%d) is greater than the maximum number of workers (%d) in your cluster profile setting. We will use %d workers instead.\n', n_jobs, max_workers, max_workers);
                fprintf("\nINFO: You may change the maximum number of workers in your cluster by running following command in the MATLAB command window:\n\n");
                fprintf("    myCluster = parcluster(); myCluster.NumWorkers = <new_number_of_workers>; saveProfile(myCluster);\n");
            end
        end
        n_jobs = min(n_jobs, max_workers);
        if ~silent
            fprintf('\nINFO: Opening a parallel pool with %d workers...\n\n', n_jobs);
            parpool(n_jobs);
        else
            evalc("parpool(n_jobs);");
        end
        success = true;
    catch ME
        fprintf('\nINFO: Failed to open a parallel pool with %d workers. Error message: %s\n', n_jobs, ME.message);
        success = false;
    end
end