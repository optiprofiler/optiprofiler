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

    % Check whether there is an existing parallel pool.
    existing_pool = gcp('nocreate');
    if ~isempty(existing_pool)
        % Borrow caller-owned resources, including smaller pools. Replacing a
        % pool can discard their worker state or interrupt unrelated work.
        if ~silent
            printOptiProfilerMessage('INFO', sprintf('Using the existing pool with up to %d workers (n_jobs=%d); the pool will not be resized or deleted.', min(n_jobs, existing_pool.NumWorkers), n_jobs));
        end
        success = true;
        return;
    end

    % Try to open a new parallel pool with maximum workers.
    try
        defaultCluster = parcluster();
        max_workers = defaultCluster.NumWorkers;
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
