function [stamp, time_stamp] = reserveExperimentDirectory(path_out, solver_names, problem_options, feature_stamp, time_stamp)
%RESERVEEXPERIMENTDIRECTORY Atomically reserve output without overwriting a run.
%   The numeric collision suffix is independent of all experiment random streams.

    if ~isfolder(path_out)
        [ok, message] = mkdir(path_out);
        if ~ok
            error('OptiProfiler:OutputDirectory', '%s', message);
        end
    end
    base_time_stamp = time_stamp;
    collision = 0;
    while true
        if collision > 0
            time_stamp = sprintf('%s_%03d', base_time_stamp, collision);
        end
        stamp = createStamp(solver_names, problem_options, feature_stamp, time_stamp, path_out);
        if exclusiveDirectory(fullfile(path_out, stamp))
            return;
        else
            collision = collision + 1;
        end
    end
end
