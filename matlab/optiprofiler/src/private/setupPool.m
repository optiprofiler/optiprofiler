function [success, reason] = setupPool(n_jobs, silent)
%SETUPPOOL set up the parallel pool by the number of jobs.
%   REASON names why no pool can be used ('' when SUCCESS is true).

    reason = '';
    % Check whether the Parallel Computing Toolbox is available.
    if exist('parpool', 'file') ~= 2
        if ~silent
            printOptiProfilerMessage('WARNING', 'Parallel Computing Toolbox is not available; the problems are solved sequentially.');
            fprintf("\n");
        end
        success = false;
        reason = 'parallel_computing_toolbox_unavailable';
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
            fprintf("\n");
        else
            evalc("parpool(n_jobs);");
        end
        success = true;
    catch ME
        if ~silent
            fprintf('\n');
            printOptiProfilerMessage('WARNING', sprintf('Failed to open a parallel pool with %d workers; the problems are solved sequentially.', n_jobs));
            printOptiProfilerMessage('WARNING', sprintf('Error message: %s', shortenMessageForLog(ME.message)));
        end
        success = false;
        reason = 'parallel_pool_start_failed';
    end
end
