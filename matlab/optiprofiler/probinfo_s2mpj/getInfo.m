function getInfo()

    currentDir = pwd;
    problem_path = fullfile(currentDir, 'matlab_problems');
    addpath(problem_path);
    filename = 'list_of_matlab_problems';
    fid = fopen(filename, 'r');
    if fid == -1
        error('Cannot open file: %s', filename);
    end
    problem_names = textscan(fid, '%s');
    fclose(fid);
    problem_names = problem_names{1};

    % Initialize the structure to store data
    problem_data = cell(length(problem_names) + 1, 9);
    problem_data{1, 1} = 'name';
    problem_data{1, 2} = 'type';
    problem_data{1, 3} = 'dim';
    problem_data{1, 4} = 'm';
    problem_data{1, 5} = 'm_linear_ub';
    problem_data{1, 6} = 'm_linear_eq';
    problem_data{1, 7} = 'm_nonlinear_ub';
    problem_data{1, 8} = 'm_nonlinear_eq';
    problem_data{1, 9} = 'loadingtime';


    for i_problem = 2:length(problem_names) + 1
        problem_name = problem_names{i_problem - 1};
        problem_name = strrep(problem_name, '.m', '');  % Remove the .m extension
        problem_data{i_problem, 1} = problem_name;

        % Try to load the problem and ignore all the errors and messages
        fprintf('Loading problem %i: %s\n', i_problem - 1, problem_name);
        tic;
        try
            [~, p] = evalc('loader(problem_name)');
            problem_data{i_problem, 9} = toc;
            fprintf('Problem %i loaded in %f seconds\n', i_problem - 1, problem_data{i_problem, 9});
        catch exception
            toc;
            problem_data{i_problem, 1} = [problem_name, ' (error loading)'];
            problem_data{i_problem, 9} = NaN;
            fprintf('Error loading problem %i: %s\n', i_problem - 1, problem_name);
        end
        try
            switch p.type
                case 'nonlinearly constrained'
                    problem_data{i_problem, 2} = 'n';
                case 'linearly constrained'
                    problem_data{i_problem, 2} = 'l';
                case 'bound-constrained'
                    problem_data{i_problem, 2} = 'b';
                case 'unconstrained'
                    problem_data{i_problem, 2} = 'u';
            end
        catch exception
            problem_data{i_problem, 2} = 'unknown';
        end
        try
            problem_data{i_problem, 3} = p.n;
        catch exception
            problem_data{i_problem, 3} = 'unknown';
        end
        try
            problem_data{i_problem, 5} = p.m_linear_ub;
        catch exception
            problem_data{i_problem, 5} = 'unknown';
        end
        try
            problem_data{i_problem, 6} = p.m_linear_eq;
        catch exception
            problem_data{i_problem, 6} = 'unknown';
        end
        try
            problem_data{i_problem, 7} = p.m_nonlinear_ub;
        catch exception
            problem_data{i_problem, 7} = 'unknown';
        end
        try
            problem_data{i_problem, 8} = p.m_nonlinear_eq;
        catch exception
            problem_data{i_problem, 8} = 'unknown';
        end
        try
            problem_data{i_problem, 4} = p.m_linear_eq + p.m_linear_ub + p.m_nonlinear_eq + p.m_nonlinear_ub;
        catch exception
            problem_data{i_problem, 4} = 'unknown';
        end
    
    end

    save('probinfo.mat', 'problem_data');

end