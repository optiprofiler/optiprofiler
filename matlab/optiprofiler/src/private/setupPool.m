function success = setupPool(n_jobs, silent)
%SETUPPOOL set up the parallel pool by the number of jobs.

    % Check whether the Parallel Computing Toolbox is available.
    if exist('parpool', 'file') ~= 2
        if ~silent
            printOptiProfilerMessage('INFO', 'Parallel Computing Toolbox is not available. We will not use parallel computing.');
        end
        success = false;
        fprintf("\n");
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
                printOptiProfilerMessage('INFO', sprintf('There is an existing parallel pool with %d workers. We will use %d workers by the option `n_jobs`.', n_workers, n_jobs));
            end
            success = true;
            fprintf("\n");
            return;
        elseif n_workers == max_workers
            % If there are maximum workers in the cluster profile, then we will use all workers.
            if ~silent
                printOptiProfilerMessage('INFO', sprintf('There is an existing parallel pool with %d workers, which is equal to the maximum number of workers in your cluster profile setting. We will use all %d workers instead of the option `n_jobs` (%d).', n_workers, n_workers, n_jobs));
                fprintf("\n");
                printOptiProfilerMessage('INFO', 'You may change the maximum number of workers in your cluster by running following command in the MATLAB command window:');
                fprintf("\n");
                fprintf("    myCluster = parcluster(); myCluster.NumWorkers = <new_number_of_workers>; saveProfile(myCluster);\n\n");
            end
            success = true;
            fprintf("\n");
            return;
        else
            % If there are fewer workers than `n_jobs`, then we will close the existing pool and try to open a new one with more workers later.
            if ~silent
                printOptiProfilerMessage('INFO', sprintf('There is an existing parallel pool with %d workers, which is fewer than the option `n_jobs` (%d). We will close the existing pool and try to open a new one with more workers.', n_workers, n_jobs));
            end
            delete(myCluster);
        end
    end

    % Try to open a new parallel pool with maximum workers.
    try
        if n_jobs > max_workers
            if ~silent
                fprintf('\n');
                printOptiProfilerMessage('INFO', sprintf('The option `n_jobs` (%d) is greater than the maximum number of workers (%d) in your cluster profile setting. We will use %d workers instead.', n_jobs, max_workers, max_workers));
                fprintf("\n");
                printOptiProfilerMessage('INFO', 'You may change the maximum number of workers in your cluster by running following command in the MATLAB command window:');
                fprintf("\n");
                fprintf("    myCluster = parcluster(); myCluster.NumWorkers = <new_number_of_workers>; saveProfile(myCluster);\n");
            end
        end
        n_jobs = min(n_jobs, max_workers);
        if ~silent
            printOptiProfilerMessage('INFO', sprintf('Starting a parallel pool with %d workers...', n_jobs));
            fprintf('\n');
            parpool(n_jobs);
        else
            evalc("parpool(n_jobs);");
        end
        success = true;
        fprintf("\n");
    catch ME
        fprintf('\n');
        printOptiProfilerMessage('INFO', sprintf('Failed to open a parallel pool with %d workers.', n_jobs));
        printOptiProfilerMessage('INFO', sprintf('Error message: %s', shortenMessageForLog(ME.message)));
        success = false;
    end
end
