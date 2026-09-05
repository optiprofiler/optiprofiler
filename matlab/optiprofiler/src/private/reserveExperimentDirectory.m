function [stamp, time_stamp] = reserveExperimentDirectory(path_out, solver_names, problem_options, feature_stamp, time_stamp)
%RESERVEEXPERIMENTDIRECTORY Atomically reserve output without overwriting a run.
%   The numeric collision suffix is independent of all experiment random streams.

    if ~usejava('jvm')
        error('OptiProfiler:OutputDirectoryJVMRequired', ...
            'Saving experiments safely requires MATLAB with the JVM enabled. Do not use -nojvm when saving benchmark output.');
    end
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
        directory = java.io.File(fullfile(path_out, stamp));
        % File.mkdir returns false when another process already created the
        % directory, unlike MATLAB mkdir's successful "already exists" result.
        if directory.mkdir()
            return;
        elseif directory.exists()
            collision = collision + 1;
        else
            error('OptiProfiler:OutputDirectory', 'Cannot create experiment directory ''%s''.', char(directory.getPath()));
        end
    end
end
