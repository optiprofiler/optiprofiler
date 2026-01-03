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
    matcutest_dir = fullfile(plib_dir, 'matcutest'); % Directory containing tools (interfaces) for MatCUTEst
    
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
    [action, ~, wrong_input] = parse_input(varargin);
    
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
        
        % Define source (Python) and destination (MATLAB) paths
        s2mpj_src_python = fullfile(setup_dir, 'python', 'optiprofiler', 'problem_libs', 's2mpj', 'src');
        s2mpj_dest_matlab = s2mpj_dir;
        s2mpj_dest_src = fullfile(s2mpj_dest_matlab, 'src');

        % Create destination source directory if it does not exist
        if ~exist(s2mpj_dest_src, 'dir')
            mkdir(s2mpj_dest_src);
        end
        
        % Attempt to copy S2MPJ source from Python directory
        if exist(s2mpj_dest_src, 'dir') && length(dir(s2mpj_dest_src)) > 2
            fprintf('S2MPJ already exists in the MATLAB directory. Skipping copy.\n');
        elseif exist(s2mpj_src_python, 'dir')
            fprintf('Copying S2MPJ from Python directory...\n');
            try
                % Copy content of s2mpj_src_python to s2mpj_dest_src
                copyfile(s2mpj_src_python, s2mpj_dest_src);
                fprintf('S2MPJ copied successfully.\n');
            catch ME
                fprintf('WARNING: Failed to copy S2MPJ automatically.\n');
                fprintf('Error: %s\n', ME.message);
                fprintf('Please manually copy the contents of:\n  %s\nto:\n  %s\n', s2mpj_src_python, s2mpj_dest_src);
            end
        else
            % If source doesn't exist, we can't copy. User might have moved things or it's a fresh clone without submodules.
            fprintf('WARNING: S2MPJ source not found at %s.\n', s2mpj_src_python);
            % Try to clone S2MPJ from its Git repository.
            fprintf('Cloning the S2MPJ repository from GitHub... This may take a few minutes.\n');
            
            clone_cmd = sprintf('git clone https://github.com/GrattonToint/S2MPJ.git "%s"', s2mpj_dest_src);
            status = system(clone_cmd);
            if status == 0
                fprintf('S2MPJ cloned successfully.\n');
            else
                fprintf('ERROR: Failed to clone S2MPJ repository. Please download S2MPJ manually from:\n');
                fprintf('https://github.com/GrattonToint/S2MPJ\n');
                fprintf("and name the directory 'src' under:\n  %s\n", s2mpj_dest_matlab);
                return; % Stop setup if cloning fails
            end
        end
        
        % Verify S2MPJ installation
        % We check if the src directory is populated (has more than just . and ..)
        s2mpj_files = dir(s2mpj_dest_src);
        if length(s2mpj_files) <= 2
            fprintf('ERROR: S2MPJ library not found in %s.\n', s2mpj_dest_src);
            fprintf('Please ensure S2MPJ is properly installed or copied manually.\n');
            return; % Stop setup if S2MPJ is missing
        end
        
        
        % =================================================================
        % 2. MatCUTEst Setup (Linux only)
        % =================================================================
        fprintf('\n--- Setting up MatCUTEst ---\n\n');
        
        paths_to_add = {src_dir, s2mpj_dir};
        
        if isunix() && ~ismac()
            % Check if MatCUTEst is already installed
            % Strategy: Check if 'matcutest' is in path or help returns something
            is_matcutest_installed = false;
            
            % Check 1: exist function
            if exist('matcutest', 'file') == 2 || exist('matcutest', 'file') == 3
                is_matcutest_installed = true;
            else
                % Check 2: help command
                try
                    help_str = help('matcutest');
                    if ~isempty(help_str)
                        is_matcutest_installed = true;
                    end
                catch
                    % Ignore errors from help
                end
            end
            
            if is_matcutest_installed
                fprintf('MatCUTEst is already installed.\n');
                paths_to_add{end+1} = matcutest_dir;
            else
                fprintf('MatCUTEst is NOT detected.\n');
                
                % Ask user if they want to install it
                user_response = input('Do you want to download and install MatCUTEst? (y/n): ', 's');
                user_response = strtrim(lower(user_response));
                
                if strcmpi(user_response, 'y')
                    % Define paths
                    % matcutest_src_dir: The final destination for the installed MatCUTEst files.
                    matcutest_src_dir = fullfile(matcutest_dir, 'src');
                    % matcutest_clone_dir: A temporary directory to clone the repository into.
                    % We use a temporary directory to avoid conflicts with the final 'src' directory
                    % and to allow for a clean installation process as recommended by MatCUTEst.
                    matcutest_clone_dir = fullfile(matcutest_dir, 'matcutest_compiled_git');
                    
                    % Clean up any previous incomplete attempts
                    if exist(matcutest_clone_dir, 'dir')
                        rmdir(matcutest_clone_dir, 's');
                    end
                    if exist(matcutest_src_dir, 'dir')
                        rmdir(matcutest_src_dir, 's');
                    end
                    
                    fprintf('Cloning MatCUTEst repository... This may take a few minutes.\n');
                    % Git clone command
                    % We clone into a temporary directory first.
                    clone_cmd = sprintf('git clone https://github.com/matcutest/matcutest_compiled.git "%s"', matcutest_clone_dir);
                    
                    status = system(clone_cmd);
                    
                    if status == 0
                        fprintf('MatCUTEst cloned successfully.\n');
                        
                        % Run internal installation
                        % The official installation command is:
                        %   rm -rf matcutest_compiled && git clone ... && matlab -batch "cd matcutest_compiled; install;" && rm -rf matcutest_compiled
                        %
                        % We adapt this by:
                        % 1. Cloning to 'matcutest_clone_dir'.
                        % 2. Running 'install(matcutest_src_dir)' to install to our desired location.
                        % 3. Removing 'matcutest_clone_dir'.
                        
                        fprintf('Running MatCUTEst installation script...\n');
                        current_dir = pwd;
                        try
                            cd(matcutest_clone_dir);
                            % Check for install.m
                            if exist('install.m', 'file')
                                % Run the install script/function.
                                % We pass the destination directory to the install function.
                                try
                                    install(matcutest_src_dir);
                                    fprintf('MatCUTEst installed successfully to %s.\n', matcutest_src_dir);
                                    
                                    % Add the directory to the path
                                    paths_to_add{end+1} = matcutest_dir;
                                    
                                    % Clean up the cloned repository
                                    cd(current_dir); % Move out before deleting
                                    fprintf('Cleaning up cloned repository...\n');
                                    rmdir(matcutest_clone_dir, 's');
                                catch ME_install
                                    fprintf('WARNING: MatCUTEst installation script failed: %s\n', ME_install.message);
                                    cd(current_dir);
                                end
                            else
                                fprintf('WARNING: install.m not found in MatCUTEst repository.\n');
                                cd(current_dir);
                            end
                        catch ME
                            fprintf('WARNING: An error occurred during MatCUTEst installation: %s\n', ME.message);
                            cd(current_dir); % Ensure we switch back even on error
                        end
                    else
                        fprintf('WARNING: Failed to clone MatCUTEst. Please install manually.\n');
                    end
                else
                    fprintf('Skip MatCUTEst installation.\n');
                end
            end
        else
            fprintf('Not a Linux system. MatCUTEst is not supported and will be skipped.\n');
        end
        
        
        % =================================================================
        % 3. Path Configuration & Persistence
        % =================================================================
        fprintf('\n--- Finalizing Setup ---\n');
        
        paths_saved = add_save_path(paths_to_add, package_name);

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


function [action, options, wrong_input] = parse_input(argin)
    %PARSE_INPUT parses the input to the setup script.
    
    action_list = {'install', 'uninstall'};
    action = 'install';
    options = []; % Initialize options
    wrong_input = false;
    
    % Start the parsing to set `input_string` and `options`.
    input_string = 'install';  % Default value for `input_string`.
    if length(argin) > 2
        fprintf('\nSetup accepts at most two inputs.\n\n');
        wrong_input = true;
    elseif length(argin) == 1
        if ischarstr(argin{1})
            input_string = argin{1};
        elseif isa(argin{1}, 'struct') || isempty(argin{1})
            options = argin{1};
        else
            fprintf('\nThe input to setup should be a string and/or a structure.\n\n');
            wrong_input = true;
        end
    elseif length(argin) == 2
        if (ischarstr(argin{1})) && (isa(argin{2}, 'struct') || isempty(argin{2}))
            input_string = argin{1};
            options = argin{2};
        elseif (ischarstr(argin{2})) && (isa(argin{1}, 'struct') || isempty(argin{1}))
            input_string = argin{2};
            options = argin{1};
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
        wrong_input = true;
    end
    
    return

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function iscs = ischarstr(x)
    %ISCHARSTR checks whether an input is a `char` or `string`
    
    iscs = (isa(x, 'char') || isa(x, 'string'));
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

    % Check whether MatCUTEst is installed inside OptiProfiler
    if isunix() && ~ismac()
        matcutest_src_dir = fullfile(matcutest_dir, 'src');
        if exist(matcutest_src_dir, 'dir')
            % TODO: Uninstall MatCUTEst properly if needed
        end
    end

    % Try removing the paths possibly added by OptiProfiler
    orig_warning_state = warning;
    warning('off', 'MATLAB:rmpath:DirNotFound'); 
    warning('off', 'MATLAB:SavePath:PathNotSaved'); 
    
    % Standard removal of known paths
    rmpath(src_dir, s2mpj_dir, matcutest_dir);
    
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
    to_be_removed = {src_dir, s2mpj_dir, matcutest_dir};
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
        
        % 2. Robust removal: Remove any lines containing the package stamp
        try
            % Read the file content
            file_content = fileread(user_startup);
            lines = strsplit(file_content, '\n');
            
            % Identify lines with the stamp
            stamp_pattern = ['% ', path_string_stamp];
            lines_to_keep = ~contains(lines, stamp_pattern);
            
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
