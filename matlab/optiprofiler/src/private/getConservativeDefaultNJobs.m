function n_jobs = getConservativeDefaultNJobs(n_workers)
%GETCONSERVATIVEDEFAULTNJOBS Get a conservative default number of workers.

    if n_workers <= 1
        n_jobs = 1;
    else
        n_jobs = max(2, floor(n_workers / 2));
    end
end
