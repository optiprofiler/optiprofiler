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
%   Setup borrows pre-existing paths and records only its own additions.
%   Uninstall without an ownership record preserves legacy paths; it never
%   guesses ownership from directory prefixes or startup-file comments.
%   An interrupted cross-file update may leave an unowned startup entry;
%   setup diagnoses and preserves that entry rather than claiming it.
%   setup(struct('install_matcutest', true))  % Set up MatCUTEst without prompting
%   setup(struct('install_solar', true))  % Download and set up the optional SOLAR MATLAB adapter
%   setup uninstall  % Uninstall the package
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   REMARKS:
%
%   To run this script, you need to have write access to the directory that
%   contains this script and its subdirectories.
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
    
    % The full paths to several directories needed for the setup.
    setup_dir = fileparts(mfilename('fullpath')); % The directory containing this setup script
    
    % Define the new directory structure
    mat_dir = fullfile(setup_dir, 'matlab'); % Matlab directory
    optiprofiler_dir = fullfile(mat_dir, 'optiprofiler'); % Directory containing the package
    src_dir = fullfile(optiprofiler_dir, 'src'); % Directory containing the source code of the package
    plib_dir = fullfile(optiprofiler_dir, 'problem_libs'); % Directory containing problem libraries
    s2mpj_dir = fullfile(plib_dir, 's2mpj'); % Directory containing S2MPJ
    % Pin S2MPJ to the commit validated with this OptiProfiler release.
    s2mpj_repo_url = 'https://github.com/optiprofiler/s2mpj_matlab.git';
    s2mpj_commit = '8d8ae2e996293b2a373db43d55a98700d16a8909';
    matcutest_dir = fullfile(plib_dir, 'matcutest'); % Directory containing tools (interfaces) for MatCUTEst
    % Pin MatCUTEst to the commit validated with this OptiProfiler release.
    matcutest_repo_url = 'https://github.com/optiprofiler/matcutest.git';
    matcutest_commit = '605d5e5b20e63d98cb7a49fc74eafd49f6a7cb80';
    solar_dir = fullfile(plib_dir, 'solar'); % Local directory for the optional SOLAR MATLAB adapter
    % Pin SOLAR to the commit validated with this OptiProfiler release.
    solar_repo_url = 'https://github.com/optiprofiler/solar_matlab.git';
    solar_commit = '7ac9439e97b91de2eb1196677bcdb858663e0bf4';
    
    % We need write access to `setup_dir` (and its subdirectories). Return if we do not have it.
    % N.B.: This checking is NOT perfect because of the following --- but it is better than nothing.
    % 1. `fileattrib` may not reflect the attributes correctly, particularly on Windows. See
    % https://www.mathworks.com/matlabcentral/answers/296657-how-can-i-check-if-i-have-read-or-write-access-to-a-directory
    % 2. Even if we have write access to `setup_dir`, we may not have the same access to its subdirectories.
    [~, attribute] = fileattrib(setup_dir);
    if ~attribute.UserWrite
        fprintf('\nSorry, we cannot continue because we do not have write access to\n\n%s\n\n', setup_dir);
        return
    end
    
    % Parse the input.
    [action, options, wrong_input] = parse_input(varargin);
    
    % Exit if wrong input detected. Error messages have been printed during the parsing.
    if wrong_input
        error('setup: The input is invalid.');
    end
    
    % Uninstall the package if requested.
    if strcmp(action, 'uninstall')
        uninstall_optiprofiler(package_name);
        return
    end
    
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
            fprintf('S2MPJ detected at %s.\n', s2mpj_dest_matlab);
            verify_existing_git_repository_commit(s2mpj_dest_matlab, s2mpj_commit, 'S2MPJ');
        else
            fprintf('S2MPJ not found. Cloning the repository...\n');
            clone_git_repository_at_commit(s2mpj_repo_url, s2mpj_dest_matlab, s2mpj_commit, 'S2MPJ');
        end
        
        
        % =================================================================
        % 2. MatCUTEst Setup (Linux only)
        % =================================================================
        fprintf('\n--- Setting up MatCUTEst ---\n\n');
        
        paths_to_add = {src_dir, s2mpj_dir};
        path_owners = {setup_dir, setup_dir};
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
            if is_matcutest_installed && is_matcutest_dir_populated
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
                    verify_existing_git_repository_commit(matcutest_dir, matcutest_commit, 'MatCUTEst');
                else
                    fprintf('Cloning MatCUTEst repository (optiprofiler fork)...\n');
                    clone_git_repository_at_commit(matcutest_repo_url, matcutest_dir, matcutest_commit, 'MatCUTEst');
                end
                
                % 5. Add the MatCUTEst directory to the path list
                paths_to_add{end+1} = matcutest_dir;
                path_owners{end+1} = matcutest_dir;
                
                % Use only the managed adapter entry point, never native setup.
                runtime_parent = fullfile(matcutest_dir, 'src');
                expected_tools = fullfile(runtime_parent, 'matcutest', 'mtools', 'src');
                user_runtime = is_matcutest_installed && ...
                    ~strcmp(fileparts(which('macup')), expected_tools);
                if user_runtime
                    fprintf('Using the existing user-managed MatCUTEst runtime; no paths are claimed.\n');
                else
                    [~, prepare_runtime] = setup_helpers();
                    receipt = prepare_runtime(matcutest_dir, runtime_parent);
                    paths_to_add = [paths_to_add, receipt.runtime_paths];
                    path_owners = [path_owners, repmat({matcutest_dir}, size(receipt.runtime_paths))];
                    fprintf('MatCUTEst managed runtime prepared at %s.\n', receipt.runtime_root);
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
            fprintf('SOLAR MATLAB adapter detected at %s. Skipping setup query.\n', solar_dir);
            verify_existing_git_repository_commit(solar_dir, solar_commit, 'SOLAR MATLAB adapter');
            proceed_with_solar = true;
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
            if exist(plib_dir, 'dir') ~= 7
                mkdir(plib_dir);
            end
            clone_git_repository_at_commit(solar_repo_url, solar_dir, solar_commit, 'SOLAR MATLAB adapter');
            is_solar_dir_populated = true;
        elseif ~proceed_with_solar
            fprintf('Skipping SOLAR MATLAB adapter setup.\n');
        end

        if is_solar_dir_populated && proceed_with_solar
            paths_to_add{end+1} = solar_dir;
            path_owners{end+1} = solar_dir;
            fprintf('SOLAR MATLAB adapter will be added to the MATLAB path.\n');
            fprintf('The adapter vendors a slim SOLAR runtime under LGPL-2.1; see its README and runtime manifest for details.\n');
        end
        
        
        % =================================================================
        % 4. Path Configuration & Persistence
        % =================================================================
        fprintf('\n--- Finalizing Setup ---\n');
        
        paths_saved = add_save_path(paths_to_add, path_owners);

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


function paths_saved = add_save_path(path_strings, owners)
    % Only explicit additions are recorded; existing caller paths are borrowed.
    [manage_paths, ~] = setup_helpers();
    context = fileparts(mfilename('fullpath'));
    paths_saved = false(numel(path_strings), 1);
    groups = unique(owners, 'stable');
    for k = 1:numel(groups)
        indices = strcmp(owners, groups{k});
        paths_saved(indices) = manage_paths('add', context, groups{k}, path_strings(indices));
    end
end

function [manage_paths, prepare_runtime] = setup_helpers()
    % Bind private helpers without adding a private directory to MATLAB path.
    original = pwd;
    cleanup = onCleanup(@() cd(original)); %#ok<NASGU>
    private_dir = fullfile(fileparts(mfilename('fullpath')), ...
        'matlab', 'optiprofiler', 'src', 'private');
    cd(private_dir);
    manage_paths = @setupPathOwnership;
    prepare_runtime = @prepareMatcutestRuntime;
end

function clone_git_repository_at_commit(repo_url, repo_dir, commit_hash, repo_name)
    %CLONE_GIT_REPOSITORY_AT_COMMIT clones and verifies a frozen problem-library checkout.

    clone_cmd = sprintf('git clone "%s" "%s"', repo_url, repo_dir);
    fprintf('Executing: %s\n', clone_cmd);
    status = system(clone_cmd);
    if status ~= 0
        remove_failed_clone_directory(repo_dir, repo_name);
        error('setup:ProblemLibraryCloneFailed', ...
            'Failed to clone %s from %s.', repo_name, repo_url);
    end

    checkout_cmd = sprintf('git -C "%s" checkout --detach %s', repo_dir, commit_hash);
    status = system(checkout_cmd);
    if status ~= 0
        remove_failed_clone_directory(repo_dir, repo_name);
        error('setup:ProblemLibraryCheckoutFailed', ...
            'Failed to checkout %s at pinned commit %s.', repo_name, commit_hash);
    end

    [verified, actual_commit] = git_repository_is_at_commit(repo_dir, commit_hash);
    if ~verified
        remove_failed_clone_directory(repo_dir, repo_name);
        error('setup:ProblemLibraryVerificationFailed', ...
            '%s checkout is at %s instead of pinned commit %s.', ...
            repo_name, actual_commit, commit_hash);
    end

    fprintf('%s cloned at pinned commit %s.\n', repo_name, commit_hash);

    return

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function verify_existing_git_repository_commit(repo_dir, commit_hash, repo_name)
    %VERIFY_EXISTING_GIT_REPOSITORY_COMMIT checks a user-provided checkout without modifying it.

    git_dir = fullfile(repo_dir, '.git');
    if exist(git_dir, 'dir') ~= 7 && exist(git_dir, 'file') ~= 2
        warning('setup:ProblemLibraryNotGitRepository', ...
            '%s directory is not a git repository. Using it as-is.', repo_name);
        return
    end

    [verified, actual_commit] = git_repository_is_at_commit(repo_dir, commit_hash);
    if verified
        fprintf('%s repository is at pinned commit %s.\n', repo_name, commit_hash);
    else
        warning('setup:ProblemLibraryCommitMismatch', ...
            '%s repository is at %s instead of pinned commit %s. Using it as-is.', ...
            repo_name, actual_commit, commit_hash);
    end

    return

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [verified, actual_commit] = git_repository_is_at_commit(repo_dir, commit_hash)
    %GIT_REPOSITORY_IS_AT_COMMIT reports whether REPO_DIR is exactly at COMMIT_HASH.

    head_cmd = sprintf('git -C "%s" rev-parse HEAD', repo_dir);
    [status, head_output] = system(head_cmd);
    actual_commit = strtrim(head_output);
    if status ~= 0 || isempty(actual_commit)
        actual_commit = '<unknown>';
        verified = false;
    else
        verified = strcmp(actual_commit, commit_hash);
    end

    return

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function remove_failed_clone_directory(repo_dir, repo_name)
    %REMOVE_FAILED_CLONE_DIRECTORY removes a directory created by a failed setup clone.

    if exist(repo_dir, 'dir') == 7
        try
            rmdir(repo_dir, 's');
        catch exception
            warning('setup:ProblemLibraryCleanupFailed', ...
                'Could not remove the incomplete %s directory %s: %s', ...
                repo_name, repo_dir, exception.message);
        end
    end

    return

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [action, options, wrong_input] = parse_input(argin)
    %PARSE_INPUT parses the input to the setup script.
    
    action_list = {'install', 'uninstall'};
    action = 'install';
    option_list = {'install_matcutest', 'install_solar'};
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

    install_only_fields = intersect(option_fields, option_list);
    if strcmp(action, 'uninstall') && ~isempty(install_only_fields)
        fprintf('\nThe setup option(s) %s only apply to `install`; they are not used with `uninstall`.\n\n', strjoin(install_only_fields, ', '));
        wrong_options = true;
    end

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


function uninstall_optiprofiler(~)
    % Keep registry and all data; remove only setup's recorded exact paths.
    context = fileparts(mfilename('fullpath'));
    [manage_paths, ~] = setup_helpers();
    found = manage_paths('remove-context', context);
    if ~found
        fprintf('No setup ownership record exists. Existing user paths and startup bytes were preserved.\n');
    else
        fprintf('Setup-owned paths removed; borrowed paths, registrations, sources and data preserved.\n');
    end
end
