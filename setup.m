function setup(varargin)
%SETUP sets the package up for MATLAB.
%
%   This file is based on
%       https://github.com/libprima/prima/blob/main/setup.m
%   which is written by Zaikun Zhang.
%
%   This script can be called in the following ways.
%
%   setup  % Add the paths needed to use the package
%   setup(struct('install_matcutest', true))  % Set up MatCUTEst without prompting
%   setup(struct('install_solar', true))  % Download and set up the optional SOLAR MATLAB adapter
%   setup(struct('problem_library_root', '/path/to/libraries'))  % Relocate optional libraries
%   setup uninstall  % Uninstall the package
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   REMARKS:
%
%   S2MPJ is bundled with OptiProfiler. Optional problem libraries are stored
%   outside the OptiProfiler source tree under the MATLAB preferences directory
%   by default. Use the `problem_library_root` option to choose another writable
%   location.
%
%   ***********************************************************************
%   Authors:
%           Cunxin HUANG (cun-xin.huang@connect.polyu.hk)
%           Tom M. RAGONNEAU (t.ragonneau@gmail.com)
%           Zaikun ZHANG (zhangzaikun@mail.sysu.edu.cn)
%           Department of Mathematics,
%           Sun Yat-sen University
%   ***********************************************************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Attribute: public (can be called directly by users)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % setup starts
    
    % Name of the package. It will be used as a stamp to be included in the path_string. Needed only
    % if `savepath` fails.
    package_name = 'optiprofiler';
    
    % Check the version of MATLAB.
    if verLessThan('matlab', '9.11')   % MATLAB R2021b = MATLAB 9.11
        fprintf('\nSorry, this package does not support MATLAB R2021a or earlier releases.\n\n');
        return
    end
    
    % The full path to the repository root.
    setup_dir = fileparts(mfilename('fullpath')); % The directory containing this setup script

    % Parse the input before constructing optional-library paths because the
    % caller may relocate them with `problem_library_root`.
    [action, options, wrong_input] = parse_input(varargin);
    if wrong_input
        error('setup: The input is invalid.');
    end
    if strcmp(action, 'uninstall')
        uninstall_optiprofiler(package_name);
        return
    end
    
    % Define the new directory structure
    mat_dir = fullfile(setup_dir, 'matlab'); % Matlab directory
    optiprofiler_dir = fullfile(mat_dir, 'optiprofiler'); % Directory containing the package
    src_dir = fullfile(optiprofiler_dir, 'src'); % Directory containing the source code of the package
    bundled_plib_dir = fullfile(optiprofiler_dir, 'problem_libs'); % Bundled problem libraries
    s2mpj_dir = fullfile(bundled_plib_dir, 's2mpj'); % Directory containing S2MPJ
    problem_library_lock = read_matlab_problem_library_lock(setup_dir);
    s2mpj_spec = get_locked_problem_library(problem_library_lock, 's2mpj');
    matcutest_spec = get_locked_problem_library(problem_library_lock, 'matcutest');
    solar_spec = get_locked_problem_library(problem_library_lock, 'solar');
    optional_plib_dir = get_optional_problem_library_root(options);
    matcutest_dir = fullfile(optional_plib_dir, char(matcutest_spec.install_directory));
    solar_dir = fullfile(optional_plib_dir, char(solar_spec.install_directory));
    matcutest_runtime_root = fullfile(optional_plib_dir, 'runtime');
    s2mpj_repo_url = locked_repository_url(s2mpj_spec);
    matcutest_repo_url = locked_repository_url(matcutest_spec);
    solar_repo_url = locked_repository_url(solar_spec);
    s2mpj_commit = char(s2mpj_spec.commit);
    matcutest_commit = char(matcutest_spec.commit);
    solar_commit = char(solar_spec.commit);
    
    % Install the package if requested.
    if strcmp(action, 'install')
        
        % =================================================================
        % 1. S2MPJ Setup
        % =================================================================
        fprintf('\n--- Setting up S2MPJ ---\n\n');
        
        % Define destination (MATLAB) path for S2MPJ
        s2mpj_dest_matlab = s2mpj_dir;

        % Check if S2MPJ directory exists and is not empty
        if exist(s2mpj_dest_matlab, 'dir') && length(dir(s2mpj_dest_matlab)) > 2
            fprintf('S2MPJ detected at %s. No further action required.\n', s2mpj_dest_matlab);
        else
            fprintf('S2MPJ not found. Cloning the repository...\n');
            
            % Clone S2MPJ from the OptiProfiler GitHub organization to s2mpj_dest_matlab
            clone_cmd = sprintf('git clone "%s" "%s"', s2mpj_repo_url, s2mpj_dest_matlab);
            fprintf('Executing: %s\n', clone_cmd);
            
            status = system(clone_cmd);
            if status == 0
                fprintf('S2MPJ cloned successfully.\n');
            else
                fprintf('ERROR: Failed to clone S2MPJ repository.\n');
                fprintf('Please manually clone or download "s2mpj_matlab" from:\n');
                fprintf('  https://github.com/%s\n', char(s2mpj_spec.repository));
                fprintf('Then rename the folder from "s2mpj_matlab" to "s2mpj" and place it at:\n');
                fprintf('  %s\n', bundled_plib_dir);
                return; % Stop setup if cloning fails
            end
        end
        checkout_git_repository_commit_if_clean( ...
            s2mpj_dir, s2mpj_repo_url, s2mpj_commit, 'S2MPJ');
        
        
        % =================================================================
        % 2. MatCUTEst Setup (Linux only)
        % =================================================================
        fprintf('\n--- Setting up MatCUTEst ---\n\n');
        
        paths_to_add = {src_dir, s2mpj_dir};
        register_matcutest = false;
        if isunix() && ~ismac()
            % Local variable to track if we should proceed with MatCUTEst actions
            proceed_with_matcutest = false;
            skip_clone_msg = false;

            % 1. Check if MatCUTEst is already installed on the system (global/path)
            is_matcutest_installed = false;
            if exist('matcutest', 'file') == 2 || exist('matcutest', 'file') == 3
                is_matcutest_installed = true;
            else
                try
                    help_str = help('matcutest');
                    if ~isempty(help_str)
                        is_matcutest_installed = true;
                    end
                catch
                end
            end
            
            % 2. Check if the local MatCUTEst directory is populated
            is_matcutest_dir_populated = exist(matcutest_dir, 'dir') && length(dir(matcutest_dir)) > 2;

            % 3. Determine if we skip asking because we are already good to go
            % (Installed AND Populated -> Good)
            if isfield(options, 'install_matcutest') && ~options.install_matcutest
                fprintf('MatCUTEst setup disabled by setup options.\n');
            elseif is_matcutest_installed && is_matcutest_dir_populated
                fprintf('MatCUTEst is installed and local repository detected. Skipping setup query.\n');
                proceed_with_matcutest = true;
                skip_clone_msg = true; 
            else
                % Ask user if they want to download/setup MatCUTEst (includes OptiProfiler plugins)
                if isfield(options, 'install_matcutest')
                    if options.install_matcutest
                        user_response = 'y';
                    else
                        user_response = 'n';
                    end
                else
                    user_response = input('Do you want to download and install/setup MatCUTEst? (y/n): ', 's');
                end
                if strcmpi(strtrim(user_response), 'y')
                    proceed_with_matcutest = true;
                end
            end
            
            if proceed_with_matcutest
                % 4. Clone MatCUTEst repository if needed
                if exist(matcutest_dir, 'dir') && length(dir(matcutest_dir)) > 2
                    if ~skip_clone_msg
                         fprintf('MatCUTEst directory at %s seems populated. Skipping clone.\n', matcutest_dir);
                    end
                    checkout_git_repository_commit_if_clean(matcutest_dir, matcutest_repo_url, matcutest_commit, 'MatCUTEst');
                else
                    fprintf('Cloning MatCUTEst repository (optiprofiler fork)...\n');
                    ensure_directory(optional_plib_dir, 'optional problem-library root');
                    clone_cmd = sprintf('git clone "%s" "%s"', matcutest_repo_url, matcutest_dir);
                    status = system(clone_cmd);
                    
                    if status == 0
                        fprintf('MatCUTEst cloned successfully.\n');
                        checkout_git_repository_commit_if_clean(matcutest_dir, matcutest_repo_url, matcutest_commit, 'MatCUTEst');
                    else
                        fprintf('WARNING: Failed to clone MatCUTEst repository.\n');
                    end
                end
                
                % 5. Add the MatCUTEst directory to the path list
                paths_to_add{end+1} = matcutest_dir;
                register_matcutest = is_populated_directory(matcutest_dir);
                
                % 6. Install if not already installed
                if is_matcutest_installed
                    if ~skip_clone_msg
                        fprintf('MatCUTEst is already installed. Skipping installation script.\n');
                    end
                else
                    fprintf('MatCUTEst not detected. Running installation script...\n');
                    current_dir = pwd;
                    try
                        matcutest_src_dir = fullfile(matcutest_dir, 'src');
                        cd(matcutest_src_dir);
                        if exist('install.m', 'file')
                            try
                                % Run install script targeting the src directory inside the repo
                                ensure_directory(matcutest_runtime_root, 'MatCUTEst runtime root');
                                install(matcutest_runtime_root);
                                fprintf('MatCUTEst installed successfully.\n');
                            catch ME_install
                                fprintf('WARNING: MatCUTEst installation script failed: %s\n', ME_install.message);
                            end
                        else
                            fprintf('WARNING: install.m not found in %s.\n', matcutest_src_dir);
                        end
                        cd(current_dir);
                    catch ME
                        fprintf('WARNING: Error during MatCUTEst setup: %s\n', ME.message);
                        cd(current_dir);
                    end
                end
            else
                fprintf('Skipping MatCUTEst setup.\n');
            end
        else
            fprintf('Not a Linux system. MatCUTEst is not supported and will be skipped.\n');
        end


        % =================================================================
        % 3. SOLAR MATLAB adapter setup (optional)
        % =================================================================
        fprintf('\n--- Setting up SOLAR MATLAB adapter ---\n\n');

        is_solar_dir_populated = is_populated_directory(solar_dir);
        proceed_with_solar = false;
        if is_solar_dir_populated
            if isfield(options, 'install_solar') && ~options.install_solar
                fprintf('SOLAR MATLAB adapter detected at %s but disabled by setup options.\n', solar_dir);
            else
                fprintf('SOLAR MATLAB adapter detected at %s. Skipping setup query.\n', solar_dir);
                proceed_with_solar = true;
            end
        else
            if isfield(options, 'install_solar')
                if options.install_solar
                    user_response = 'y';
                else
                    user_response = 'n';
                end
            else
                user_response = input('Do you want to download/setup the SOLAR MATLAB adapter? (y/n): ', 's');
            end
            if strcmpi(strtrim(user_response), 'y')
                proceed_with_solar = true;
            end
        end

        if proceed_with_solar && ~is_solar_dir_populated
            fprintf('Cloning SOLAR MATLAB adapter (optiprofiler fork)...\n');
            ensure_directory(optional_plib_dir, 'optional problem-library root');
            clone_cmd = sprintf('git clone "%s" "%s"', solar_repo_url, solar_dir);
            status = system(clone_cmd);
            if status == 0
                fprintf('SOLAR MATLAB adapter cloned successfully.\n');
                is_solar_dir_populated = true;
            else
                fprintf('WARNING: Failed to clone SOLAR MATLAB adapter repository.\n');
                fprintf('You can manually clone or download "solar_matlab" from:\n');
                fprintf('  https://github.com/%s\n', char(solar_spec.repository));
                fprintf('Then rename the folder from "solar_matlab" to "solar" and place it at:\n');
                fprintf('  %s\n', optional_plib_dir);
            end
        elseif ~proceed_with_solar
            fprintf('Skipping SOLAR MATLAB adapter setup.\n');
        end

        if is_solar_dir_populated && proceed_with_solar
            checkout_git_repository_commit_if_clean( ...
                solar_dir, solar_repo_url, solar_commit, 'SOLAR MATLAB adapter');
            paths_to_add{end+1} = solar_dir;
            fprintf('SOLAR MATLAB adapter will be added to the MATLAB path.\n');
            fprintf('The adapter vendors a slim SOLAR runtime; see its README, license files, and runtime manifest for details.\n');
        end
        
        
        % =================================================================
        % 4. Path Configuration & Persistence
        % =================================================================
        fprintf('\n--- Finalizing Setup ---\n');
        
        paths_saved = add_save_path(paths_to_add, package_name);

        % Persist explicit optional-library registrations after the source
        % directory is available on the MATLAB path. Bundled S2MPJ is
        % resolved directly from its locked gitlink.
        if register_matcutest
            register_locked_problem_library(matcutest_spec, matcutest_dir);
        end
        if is_solar_dir_populated && proceed_with_solar
            register_locked_problem_library(solar_spec, solar_dir);
        end

        if all(paths_saved)
            fprintf('\nThe package is ready to use.\n');
            fprintf('\nYou may now try ''help optiprofiler'' for information on the usage of the package.\n');
            fprintf('\nYou may also run ''testOptiProfiler'' to test the package.\n');
            fprintf('\nA few examples showing how to use the package are provided in the directory:\n\n');
            ex_dir = fullfile(mat_dir, 'examples');
            fprintf('    %s\n\n', ex_dir);
        else
            fprintf('\n***** To use the package in other MATLAB sessions, append the following lines to your startup script. *****\n');
            fprintf('\n  (see https://www.mathworks.com/help/matlab/ref/startup.html for information):\n');
            for i = 1:length(paths_to_add)
                fprintf('    addpath(''%s'');\n', paths_to_add{i});
            end
        end

        fprintf('\n');
        return

    end
    
    % setup ends
    return

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function lock = read_matlab_problem_library_lock(setup_dir)
    %READ_MATLAB_PROBLEM_LIBRARY_LOCK reads the setup integration manifest.

    lock_path = fullfile(setup_dir, 'matlab_problem_libraries.lock');
    if exist(lock_path, 'file') ~= 2
        error('setup:ProblemLibraryLockNotFound', ...
            'Cannot find the MATLAB problem-library lock at %s.', lock_path);
    end
    try
        lock = jsondecode(fileread(lock_path));
    catch ME
        error('setup:ProblemLibraryLockInvalid', ...
            'Cannot read the MATLAB problem-library lock: %s', ME.message);
    end
    if ~isstruct(lock) || ~isfield(lock, 'problem_library_api_version') || ...
            lock.problem_library_api_version ~= 1 || ~isfield(lock, 'libraries')
        error('setup:ProblemLibraryLockInvalid', ...
            'The MATLAB problem-library lock must use API version 1.');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function entry = get_locked_problem_library(lock, name)
    %GET_LOCKED_PROBLEM_LIBRARY returns one named lock entry.

    entries = lock.libraries;
    entry = [];
    for i_entry = 1:numel(entries)
        if iscell(entries)
            candidate = entries{i_entry};
        else
            candidate = entries(i_entry);
        end
        if strcmp(char(candidate.name), name)
            entry = candidate;
            break
        end
    end
    if isempty(entry)
        error('setup:ProblemLibraryLockInvalid', ...
            'The MATLAB problem-library lock has no entry for %s.', name);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function repo_url = locked_repository_url(entry)
    %LOCKED_REPOSITORY_URL returns the HTTPS clone URL for one lock entry.

    repo_url = sprintf('https://github.com/%s.git', char(entry.repository));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function register_locked_problem_library(entry, root)
    %REGISTER_LOCKED_PROBLEM_LIBRARY persists one optional setup result.

    collect_info_function = entry.collect_info_function;
    if isempty(collect_info_function)
        collect_info_function = '';
    else
        collect_info_function = char(collect_info_function);
    end
    registration = struct( ...
        'name', char(entry.name), ...
        'api_version', double(entry.api_version), ...
        'role', char(entry.role), ...
        'root', root, ...
        'select_function', char(entry.select_function), ...
        'load_function', char(entry.load_function), ...
        'collect_info_function', collect_info_function, ...
        'check_available_function', '', ...
        'platforms', {cellstr(entry.platforms)});
    registerProblemLibrary(registration);
    fprintf('Registered problem library "%s" at %s.\n', ...
        registration.name, registration.root);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function paths_saved = add_save_path(path_strings, path_string_stamp)
    %ADD_SAVE_PATH adds the paths indicated by PATH_STRINGS to the MATLAB path and then tries saving the paths.
    % PATH_STRING_STAMP is a stamp used when writing PATH_STRINGS to the user's startup.m file, which is
    % needed only if `savepath` fails.
    
    paths_saved = false(length(path_strings), 1);
    for i_path = 1:length(path_strings)
        path_string = path_strings{i_path};

        if ~exist(path_string, 'dir')
            warning('The directory %s does not exist. Skipping.', path_string);
            continue;
        end
        addpath(path_string);
    
        % Try saving the path in the system path-defining file at sys_pathdef.
        orig_warning_state = warning;
        warning('off', 'MATLAB:SavePath:PathNotSaved'); 
        sys_pathdef = fullfile(matlabroot(), 'toolbox', 'local', 'pathdef.m');
        paths_saved(i_path) = (savepath(sys_pathdef) == 0);
        warning(orig_warning_state); 
        
        % If path not saved, try editing the startup.m of this user.
        if ~paths_saved(i_path) && numel(userpath) > 0
            user_startup = fullfile(userpath, 'startup.m');
            add_path_string = sprintf('addpath(''%s'');', path_string);
            full_add_path_string = sprintf('%s\t%s %s', add_path_string, '%', path_string_stamp);
        
            % Check if already exists
            if exist(user_startup, 'file')
                startup_text_cells = regexp(fileread(user_startup), '\n', 'split');
                paths_saved(i_path) = any(strcmp(startup_text_cells, full_add_path_string));
            end
        
            if ~paths_saved(i_path)
                % Check for empty last line
                if exist(user_startup, 'file')
                    startup_text_cells = regexp(fileread(user_startup), '\n', 'split');
                    last_line_empty = isempty(startup_text_cells) || (isempty(startup_text_cells{end}) && ...
                        isempty(startup_text_cells{max(1, end-1)}));
                else
                    last_line_empty = true;
                end
        
                file_id = fopen(user_startup, 'a');
                if file_id ~= -1 
                    if ~last_line_empty
                        fprintf(file_id, '\n');
                    end
                    fprintf(file_id, '%s', full_add_path_string);
                    fclose(file_id);
                    % Verify
                    if exist(user_startup, 'file')
                        startup_text_cells = regexp(fileread(user_startup), '\n', 'split');
                        paths_saved(i_path) = any(strcmp(startup_text_cells, full_add_path_string));
                    end
                end
            end
        end
        
    end
    
    return

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function checkout_git_repository_commit_if_clean(repo_dir, repo_url, commit_hash, repo_name)
    %CHECKOUT_GIT_REPOSITORY_COMMIT_IF_CLEAN pins an optional problem library checkout.

    git_dir = fullfile(repo_dir, '.git');
    if exist(git_dir, 'dir') ~= 7 && exist(git_dir, 'file') ~= 2
        fprintf('%s directory is not a git repository. Using it as-is.\n', repo_name);
        return
    end

    head_cmd = sprintf('git -C "%s" rev-parse HEAD', repo_dir);
    [status, head_output] = system(head_cmd);
    if status == 0 && strcmp(strtrim(head_output), commit_hash)
        fprintf('%s repository is already at pinned commit %s.\n', repo_name, commit_hash);
        return
    end

    status_cmd = sprintf('git -C "%s" status --porcelain', repo_dir);
    [status, status_output] = system(status_cmd);
    if status ~= 0
        fprintf('WARNING: Could not inspect %s repository status. Using it as-is.\n', repo_name);
        return
    end

    if ~isempty(strtrim(status_output))
        fprintf('%s repository has local changes. Skipping checkout of pinned commit %s.\n', repo_name, commit_hash);
        return
    end

    fetch_cmd = sprintf('git -C "%s" fetch --tags "%s"', repo_dir, repo_url);
    status = system(fetch_cmd);
    if status ~= 0
        fprintf('WARNING: Could not fetch %s repository. Using it as-is.\n', repo_name);
        return
    end

    checkout_cmd = sprintf('git -C "%s" checkout --detach %s', repo_dir, commit_hash);
    status = system(checkout_cmd);
    if status ~= 0
        fprintf('WARNING: Could not checkout %s pinned commit %s. Using repository as-is.\n', repo_name, commit_hash);
        return
    end

    clean_cmd = sprintf('git -C "%s" reset --hard %s', repo_dir, commit_hash);
    status = system(clean_cmd);
    if status ~= 0
        fprintf('WARNING: Could not reset %s repository to pinned commit %s. Using repository as-is.\n', repo_name, commit_hash);
        return
    end

    fprintf('%s repository checked out at pinned commit %s.\n', repo_name, commit_hash);

    return

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [action, options, wrong_input] = parse_input(argin)
    %PARSE_INPUT parses the input to the setup script.
    
    action_list = {'install', 'uninstall'};
    action = 'install';
    option_list = {'install_matcutest', 'install_solar', 'problem_library_root'};
    options = struct(); % Initialize options
    wrong_input = false;
    
    % Start the parsing to set `input_string` and `options`.
    input_string = 'install';  % Default value for `input_string`.
    if length(argin) > 2
        fprintf('\nSetup accepts at most two inputs.\n\n');
        wrong_input = true;
    elseif length(argin) == 1
        if ischarstr(argin{1})
            input_string = argin{1};
        elseif isa(argin{1}, 'struct')
            options = argin{1};
        elseif isempty(argin{1})
            options = struct();
        else
            fprintf('\nThe input to setup should be a string and/or a structure.\n\n');
            wrong_input = true;
        end
    elseif length(argin) == 2
        if (ischarstr(argin{1})) && (isa(argin{2}, 'struct') || isempty(argin{2}))
            input_string = argin{1};
            if isa(argin{2}, 'struct')
                options = argin{2};
            else
                options = struct();
            end
        elseif (ischarstr(argin{2})) && (isa(argin{1}, 'struct') || isempty(argin{1}))
            input_string = argin{2};
            if isa(argin{1}, 'struct')
                options = argin{1};
            else
                options = struct();
            end
        else
            fprintf('\nThe input to setup should be a string and/or a structure.\n\n');
            wrong_input = true;
        end
    end
    
    % Cast input_string to a character array in case it is a MATLAB string.
    input_string = lower(char(input_string));
    
    % Parse `input_string` to set `action`.
    if ismember(input_string, action_list)
        action = input_string;
    else
        fprintf('\nUnknown setup action `%s`. Valid actions are: %s.\n\n', input_string, strjoin(action_list, ', '));
        wrong_input = true;
    end

    % Validate `options` explicitly. Silently ignoring misspelled options is
    % confusing for users because setup may continue while their request has
    % no effect.
    if ~wrong_input
        [options, wrong_options] = validate_options(options, option_list, action);
        wrong_input = wrong_input || wrong_options;
    end
    
    return

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [options, wrong_options] = validate_options(options, option_list, action)
    %VALIDATE_OPTIONS validates the option structure passed to setup.

    wrong_options = false;

    if isempty(options)
        options = struct();
        return
    end

    if ~isa(options, 'struct') || ~isscalar(options)
        fprintf('\nThe options input to setup should be a scalar structure.\n\n');
        wrong_options = true;
        return
    end

    option_fields = fieldnames(options);
    unknown_fields = setdiff(option_fields, option_list);
    if ~isempty(unknown_fields)
        fprintf('\nUnknown setup option(s): %s.\n', strjoin(unknown_fields, ', '));
        fprintf('Supported setup option(s): %s.\n\n', strjoin(option_list, ', '));
        wrong_options = true;
    end

    if isfield(options, 'install_matcutest')
        if ~(islogical(options.install_matcutest) && isscalar(options.install_matcutest))
            fprintf('\nThe setup option `install_matcutest` must be a scalar logical value: true or false.\n\n');
            wrong_options = true;
        end
    end

    if isfield(options, 'install_solar')
        if ~(islogical(options.install_solar) && isscalar(options.install_solar))
            fprintf('\nThe setup option `install_solar` must be a scalar logical value: true or false.\n\n');
            wrong_options = true;
        end
    end

    if isfield(options, 'problem_library_root')
        root = options.problem_library_root;
        if ~ischarstr(root) || ~isscalar(string(root)) || isempty(strtrim(char(root)))
            fprintf('\nThe setup option `problem_library_root` must be a nonempty char or string scalar.\n\n');
            wrong_options = true;
        end
    end

    install_only_fields = intersect(option_fields, option_list);
    if strcmp(action, 'uninstall') && ~isempty(install_only_fields)
        fprintf('\nThe setup option(s) %s only apply to `install`; they are not used with `uninstall`.\n\n', strjoin(install_only_fields, ', '));
        wrong_options = true;
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function root = get_optional_problem_library_root(options)
    %GET_OPTIONAL_PROBLEM_LIBRARY_ROOT returns an absolute external cache root.

    if isfield(options, 'problem_library_root')
        root = strtrim(char(options.problem_library_root));
    else
        root = fullfile(prefdir, 'optiprofiler', 'problem_libraries');
    end
    if ~is_absolute_path(root)
        root = fullfile(pwd, root);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function tf = is_absolute_path(path_string)
    %IS_ABSOLUTE_PATH checks Unix, drive-letter, and UNC absolute paths.

    if ispc
        tf = ~isempty(regexp(path_string, '^[A-Za-z]:[\\/]', 'once')) || ...
            startsWith(path_string, '\\');
    else
        tf = startsWith(path_string, filesep);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function ensure_directory(directory, description)
    %ENSURE_DIRECTORY creates a required writable setup directory.

    if exist(directory, 'dir') ~= 7
        [status, message] = mkdir(directory);
        if ~status
            error('setup:ProblemLibraryDirectory', ...
                'Cannot create the %s at %s: %s', description, directory, message);
        end
    end
    probe = [tempname(directory), '.tmp'];
    file_id = fopen(probe, 'w');
    if file_id == -1
        error('setup:ProblemLibraryDirectory', ...
            'The %s is not writable: %s', description, directory);
    end
    fclose(file_id);
    delete(probe);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function iscs = ischarstr(x)
    %ISCHARSTR checks whether an input is a `char` or `string`
    
    iscs = (isa(x, 'char') || isa(x, 'string'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function tf = is_populated_directory(directory)
    %IS_POPULATED_DIRECTORY checks whether a directory exists and is not empty.

    tf = exist(directory, 'dir') == 7 && length(dir(directory)) > 2;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function uninstall_optiprofiler(path_string_stamp)
    %UNINSTALL_OPTIPROFILER uninstalls OptiProfiler.
    
    fprintf('\nUninstalling OptiProfiler (if it is installed) ... ');
    
    % The full path of several directories.   
    mfiledir = fileparts(mfilename('fullpath'));  % The directory where this .m file resides
    mat_dir = fullfile(mfiledir, 'matlab'); 
    optiprofiler_dir = fullfile(mat_dir, 'optiprofiler'); 
    src_dir = fullfile(optiprofiler_dir, 'src'); 
    plib_dir = fullfile(optiprofiler_dir, 'problem_libs'); 
    s2mpj_dir = fullfile(plib_dir, 's2mpj'); 
    matcutest_dir = fullfile(plib_dir, 'matcutest');
    solar_dir = fullfile(plib_dir, 'solar');
    registered_roots = get_registered_problem_library_roots();

    % We do not need to specifically uninstall S2MPJ or MatCUTEst, since they
    % only add paths to MATLAB and do not modify MATLAB installation files.
    % % Check whether MatCUTEst is installed inside OptiProfiler
    % if isunix() && ~ismac()
    %     matcutest_src_dir = fullfile(matcutest_dir, 'src');
    %     if exist(matcutest_src_dir, 'dir')
    %         % TODO: Uninstall MatCUTEst properly if needed
    %     end
    % end

    % Try removing the paths possibly added by OptiProfiler
    orig_warning_state = warning;
    warning('off', 'MATLAB:rmpath:DirNotFound'); 
    warning('off', 'MATLAB:SavePath:PathNotSaved'); 
    
    % Standard removal of known paths
    rmpath(src_dir, s2mpj_dir, matcutest_dir, solar_dir, registered_roots{:});
    
    % Robust removal: Remove any path that is a subdirectory of the package root
    % This handles cases where the folder structure might have changed or extra paths were added.
    all_paths = strsplit(path, pathsep);
    paths_to_remove_robust = {};
    for i = 1:length(all_paths)
        if startsWith(all_paths{i}, mfiledir)
            paths_to_remove_robust{end+1} = all_paths{i};
        end
    end
    if ~isempty(paths_to_remove_robust)
        rmpath(paths_to_remove_robust{:});
    end
    
    savepath;
    warning(orig_warning_state); 
    
    % Removing the line possibly added to the user startup script
    to_be_removed = [{src_dir, s2mpj_dir, matcutest_dir, solar_dir}, registered_roots];
    user_startup = fullfile(userpath,'startup.m');
    if exist(user_startup, 'file')
        % 1. Try removing specific known lines (legacy support)
        for i_path = 1:length(to_be_removed)
            add_path_string = sprintf('addpath(''%s'');', to_be_removed{i_path});
            full_add_path_string = sprintf('%s\t%s %s', add_path_string, '%', path_string_stamp);
            try
                del_str_ln(user_startup, full_add_path_string);
            catch
                % Do nothing.
            end
        end
        
        % 2. Robust removal: Remove any lines containing the package stamp OR the package root directory
        try
            % Read the file content
            file_content = fileread(user_startup);
            lines = strsplit(file_content, '\n');
            
            % Identify lines with the stamp
            stamp_pattern = ['% ', path_string_stamp];

            % Identify lines containing the root directory (handles cases where the stamp might differ or be missing)
            % This ensures we remove paths related to this package even if they were added with a different stamp (e.g., % matcutest).
            lines_to_keep = ~contains(lines, stamp_pattern) & ~contains(lines, mfiledir);
            
            if ~all(lines_to_keep)
                % Reconstruct the file content
                new_content = strjoin(lines(lines_to_keep), '\n');
                
                % Write back to file
                fid = fopen(user_startup, 'w');
                if fid ~= -1
                    fprintf(fid, '%s', new_content);
                    fclose(fid);
                end
            end
        catch
            % Ignore errors during robust cleanup
        end
    end
    
    fprintf('Done.\nYou may now remove\n\n    %s\n\nif it contains nothing you want to keep.\n\n', mfiledir);
    
    return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function roots = get_registered_problem_library_roots()
    %GET_REGISTERED_PROBLEM_LIBRARY_ROOTS reads paths persisted by the registry.

    registry_file = getenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY');
    if isempty(registry_file)
        registry_file = fullfile(prefdir, 'optiprofiler', 'problem_libraries.mat');
    end
    roots = {};
    if exist(registry_file, 'file') ~= 2
        return
    end
    try
        data = load(registry_file, 'registrations');
        if isfield(data, 'registrations') && isstruct(data.registrations)
            roots = {data.registrations.root};
            roots = roots(~cellfun(@isempty, roots));
        end
    catch
        % A stale registry should not prevent core path cleanup.
        roots = {};
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function del_str_ln(filename, string)
    %DEL_STR_LN deletes from filename all the lines that are identical to string
    
    fid = fopen(filename, 'r');  % Open file for reading.
    if fid == -1
        error('Cannot open file %s.', filename);
    end
    
    % Read the file into a cell of strings
    data = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');
    fclose(fid);
    cstr = data{1};
    
    % Remove the rows containing string
    cstr(strcmp(cstr, string)) = [];
    
    % Save the file again
    fid = fopen(filename, 'w');  % Open/create new file for writing. Discard existing contents, if any.
    if fid == -1
        error('Cannot open file %s.', filename);
    end
    fprintf(fid, '%s\n', cstr{:});
    fclose(fid);
    
    return

end
