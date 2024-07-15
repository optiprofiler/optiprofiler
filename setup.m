function setup(varargin)
    %SETUP sets the package up for MATLAB.
    %
    %   This file is based on https://github.com/libprima/prima/blob/main/setup.m, which is written by Zaikun Zhang.
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
    %   Authors:    Cunxin HUANG (cun-xin.huang@connect.polyu.hk)
    %               Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
    %               Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
    %               Department of Applied Mathematics,
    %               The Hong Kong Polytechnic University.
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
    matd = fullfile(setup_dir, 'matlab'); % Matlab directory
    optprofiler_dir = fullfile(matd, 'optiprofiler'); % Directory containing the package
    src_dir = fullfile(optprofiler_dir, 'src'); % Directory containing the source code of the package
    probinfo_dir = fullfile(optprofiler_dir, 'probinfo_s2mpj'); % Directory containing the problem information of the test set S2MPJ
    s2mpj_dir = fullfile(probinfo_dir, 'matlab_problems'); % Directory containing the test set S2MPJ
    tests_dir = fullfile(optprofiler_dir, 'tests'); % Directory containing some tests
    
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
    [action, wrong_input] = parse_input(varargin);
    
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

        paths_saved = add_save_path({src_dir, probinfo_dir, s2mpj_dir, tests_dir}, package_name);

        if paths_saved(1)
            fprintf('\nThe package is ready to use.\n');
            fprintf('\nYou may now try ''help runBenchmark'' for information on the usage of the package.\n');
            if paths_saved(2)
                fprintf('\nYou may also run ''testOptiProfiler'' to test the package on a few examples.\n');
            end
        else
            add_path_string_1 = sprintf('addpath(''%s'');', src_dir);
            add_path_string_2 = sprintf('addpath(''%s'');', probinfo_dir);
            add_path_string_3 = sprintf('addpath(''%s'');', s2mpj_dir);
            fprintf('\n***** To use the package in other MATLAB sessions, do ONE of the following. *****\n');
            fprintf('\n- Append the following line to your startup script');
            fprintf('\n  (see https://www.mathworks.com/help/matlab/ref/startup.html for information):\n');
            fprintf('\n    %s\n', add_path_string_1);
            fprintf('\n    %s\n', add_path_string_2);
            fprintf('\n    %s\n', add_path_string_3);
            fprintf('\n- OR come to the current directory and run ''setup path'' when you need the package.\n');
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
            error('The directory %s does not exist.', path_string);
        end
        addpath(path_string);
    
        % Try saving the path in the system path-defining file at sys_pathdef. If the user does not have
        % writing permission for this file, then the path will not saved.
        % N.B. Do not save the path to the pathdef.m file under userpath. This file is not loaded by default
        % at startup. See
        % https://www.mathworks.com/matlabcentral/answers/269482-is-userdata-pathdef-m-for-local-path-additions-supported-on-linux
        orig_warning_state = warning;
        warning('off', 'MATLAB:SavePath:PathNotSaved'); % Maybe we do not have the permission to save path
        sys_pathdef = fullfile(matlabroot(), 'toolbox', 'local', 'pathdef.m');
        paths_saved(i_path) = (savepath(sys_pathdef) == 0);
        warning(orig_warning_state); % Restore the behavior of displaying warnings
        
        % If path not saved, try editing the startup.m of this user. Do this only if userpath is nonempty.
        % Otherwise, we will only get a startup.m in the current directory, which will not be executed
        % when MATLAB starts from other directories. On Linux, the default value of userpath is
        % ~/Documents/MATLAB, but it will be '' if this directory does not exist. We refrain from creating
        % this directory in that case.
        if ~paths_saved(i_path) && numel(userpath) > 0
            user_startup = fullfile(userpath, 'startup.m');
            add_path_string = sprintf('addpath(''%s'');', path_string);
            full_add_path_string = sprintf('%s\t%s %s', add_path_string, '%', path_string_stamp);
        
            % First, check whether full_add_path_string already exists in user_startup or not.
            if exist(user_startup, 'file')
                startup_text_cells = regexp(fileread(user_startup), '\n', 'split');
                paths_saved(i_path) = any(strcmp(startup_text_cells, full_add_path_string));
            end
        
            if ~paths_saved(i_path)
                % We first check whether the last line of the user startup script is an empty line (or the
                % file is empty or even does not exist at all). If yes, we do not need to put a line break
                % before the path adding string.
                if exist(user_startup, 'file')
                    startup_text_cells = regexp(fileread(user_startup), '\n', 'split');
                    last_line_empty = isempty(startup_text_cells) || (isempty(startup_text_cells{end}) && ...
                        isempty(startup_text_cells{max(1, end-1)}));
                else
                    last_line_empty = true;
                end
        
                file_id = fopen(user_startup, 'a');  % Open/create file for writing. Append data to the end.
                if file_id ~= -1 % If FOPEN cannot open the file, it returns -1; We keep silent if it fails.
                    if ~last_line_empty  % The last line of user_startup is not empty
                        fprintf(file_id, '\n');  % Add a new empty line
                    end
                    fprintf(file_id, '%s', full_add_path_string);
                    fclose(file_id);
                    % Check that full_add_path_string is indeed added to user_startup.
                    if exist(user_startup, 'file')
                        startup_text_cells = regexp(fileread(user_startup), '\n', 'split');
                        paths_saved(i_path) = any(strcmp(startup_text_cells, full_add_path_string));
                    end
                end
            end
        end
        
    end
    
    % add_save_path ends
    return

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [action, wrong_input] = parse_input(argin)
    %PARSE_INPUT parses the input to the setup script.
    
    action_list = {'install', 'uninstall'};
    action = 'install';
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
    %UNINSTALL_OPTIPROFIER uninstalls OptiProfiler.
    
    fprintf('\nUninstalling OptiProfiler (if it is installed) ... ');
    
    % The full path of several directories.   
    mfiledir = fileparts(mfilename('fullpath'));  % The directory where this .m file resides
    matd = fullfile(mfiledir, 'matlab'); % Matlab directory
    optprofiler_dir = fullfile(matd, 'optiprofiler'); % Directory containing the package
    src_dir = fullfile(optprofiler_dir, 'src'); % Directory containing the source code of the package
    probinfo_dir = fullfile(optprofiler_dir, 'probinfo_s2mpj'); % Directory containing the problem information of the test set S2MPJ
    s2mpj_dir = fullfile(probinfo_dir, 'matlab_problems'); % Directory containing the test set S2MPJ
    tests_dir = fullfile(optprofiler_dir, 'tests'); % Directory containing some tests

    % Try removing the paths possibly added by OptiProfiler
    orig_warning_state = warning;
    warning('off', 'MATLAB:rmpath:DirNotFound'); % Maybe the paths were not added. We do not want to see this warning.
    warning('off', 'MATLAB:SavePath:PathNotSaved'); % Maybe we do not have the permission to save path.
    rmpath(src_dir, probinfo_dir, s2mpj_dir, tests_dir);
    savepath;
    warning(orig_warning_state); % Restore the behavior of displaying warnings
    
    % Removing the line possibly added to the user startup script
    to_be_removed = {src_dir, probinfo_dir, s2mpj_dir, tests_dir};
    user_startup = fullfile(userpath,'startup.m');
    if exist(user_startup, 'file')
        for i_path = 1:length(to_be_removed)
            add_path_string = sprintf('addpath(''%s'');', to_be_removed{i_path});
            full_add_path_string = sprintf('%s\t%s %s', add_path_string, '%', path_string_stamp);
            try
                del_str_ln(user_startup, full_add_path_string);
            catch
                % Do nothing.
            end
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