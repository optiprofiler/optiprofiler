classdef TestSetupProblemLibraryPins < matlab.unittest.TestCase
% setup brings a problem library to its locked commit or says that it cannot.
% Everything happens in a sandbox: a copy of the package with its own lock, local
% fixture repositories that stand in for GitHub, a private home and an isolated
% registry, pathdef and startup file. There is no network access, and the real
% MATLAB path and startup files are not touched.
    methods (TestClassSetup)
        function requireGit(testCase)
            testCase.assumeTrue(isunix);
            [status, ~] = system('git --version');
            testCase.assumeEqual(status, 0);
        end
    end

    methods (Test)
        function branchOnlyPinIsFetchedAndCheckedOut(testCase)
            state = sandbox(testCase);
            base = make_remote(state, 'solar_matlab', solar_files());
            clone(state, 'solar_matlab', 'solar');
            pin = branch_commit(state, 'solar_matlab');
            testCase.assertNotEqual(pin, base);
            lock_commits(state, struct('solar', pin));
            run_setup(state, struct('install_solar', true));
            testCase.verifyEqual(head(state, 'solar'), pin);
            testCase.verifyTrue(is_registered(state, 'solar'));
            testCase.verifyTrue(on_path(fullfile(state.libraries, 'solar')));
        end

        function pinHeldLocallyNeedsNoRemote(testCase)
            state = sandbox(testCase);
            make_remote(state, 'solar_matlab', solar_files());
            pin = branch_commit(state, 'solar_matlab');
            clone(state, 'solar_matlab', 'solar');
            rmdir(fullfile(state.remotes, 'optiprofiler', 'solar_matlab.git'), 's');
            lock_commits(state, struct('solar', pin));
            output = run_setup(state, struct('install_solar', true));
            testCase.verifySubstring(output, 'checked out at pinned commit');
            testCase.verifyEqual(head(state, 'solar'), pin);
        end

        function repeatedSetupStaysAtThePin(testCase)
            state = sandbox(testCase);
            pin = make_remote(state, 'solar_matlab', solar_files());
            clone(state, 'solar_matlab', 'solar');
            lock_commits(state, struct('solar', pin));
            run_setup(state, struct('install_solar', true));
            output = run_setup(state, struct('install_solar', true));
            testCase.verifySubstring(output, 'already at pinned commit');
            testCase.verifyEqual(head(state, 'solar'), pin);
        end

        function unreachablePinOfTheBundledLibraryIsAnError(testCase)
            state = sandbox(testCase);
            make_remote(state, 's2mpj_matlab', {'README.md', 'fixture'; ...
                fullfile('src', 'matlab_problems', 'keep.txt'), 'fixture'});
            s2 = fullfile(state.package, 'matlab', 'optiprofiler', 'problem_libs', 's2mpj');
            rmdir(s2, 's');
            git(state.root, sprintf('clone -q "%s" "%s"', remote_path(state, 's2mpj_matlab'), s2));
            before = rev_parse(s2);
            lock_commits(state, struct('s2mpj', unknown_commit()));
            testCase.verifyError(@() run_setup(state, struct()), 'setup:PinnedCommitUnavailable');
            testCase.verifyEqual(rev_parse(s2), before);
        end

        function unreachablePinOfAnOptionalLibraryIsAnError(testCase)
            state = sandbox(testCase);
            base = make_remote(state, 'solar_matlab', solar_files());
            clone(state, 'solar_matlab', 'solar');
            lock_commits(state, struct('solar', unknown_commit()));
            testCase.verifyError(@() run_setup(state, struct('install_solar', true)), 'setup:PinnedCommitUnavailable');
            testCase.verifyEqual(head(state, 'solar'), base);
            testCase.verifyFalse(is_registered(state, 'solar'));
        end

        function localChangesAreLeftAlone(testCase)
            state = sandbox(testCase);
            base = make_remote(state, 'solar_matlab', solar_files());
            clone(state, 'solar_matlab', 'solar');
            pin = branch_commit(state, 'solar_matlab');
            edited = fullfile(state.libraries, 'solar', 'solar_load.m');
            write_text(edited, sprintf('function p = solar_load(~)\n%% my local edit\np = [];\nend\n'));
            lock_commits(state, struct('solar', pin));
            output = run_setup(state, struct('install_solar', true));
            testCase.verifySubstring(output, 'has local changes');
            testCase.verifyEqual(head(state, 'solar'), base);
            testCase.verifySubstring(fileread(edited), 'my local edit');
            testCase.verifyTrue(is_registered(state, 'solar'));
        end

        function directoryThatIsNotARepositoryIsUsedAsIs(testCase)
            state = sandbox(testCase);
            root = fullfile(state.libraries, 'solar');
            mkdir(root);
            files = solar_files();
            for k = 1:size(files, 1)
                write_text(fullfile(root, files{k, 1}), files{k, 2});
            end
            lock_commits(state, struct('solar', unknown_commit()));
            output = run_setup(state, struct('install_solar', true));
            testCase.verifySubstring(output, 'not a git repository');
            testCase.verifyTrue(is_registered(state, 'solar'));
        end
    end
end


function state = sandbox(testCase)
    state.root = tempname;
    mkdir(state.root);
    [~, attributes] = fileattrib(state.root);
    state.root = attributes.Name;
    state.home = fullfile(state.root, 'home');
    state.remotes = fullfile(state.root, 'remotes');
    state.libraries = fullfile(state.root, 'libraries');
    state.package = fullfile(state.root, 'package');
    state.registry = fullfile(state.root, 'registry.mat');
    mkdir(state.home);
    mkdir(fullfile(state.remotes, 'optiprofiler'));
    mkdir(state.libraries);

    old_path = path;
    old_dir = pwd;
    variables = {'HOME', 'XDG_CONFIG_HOME', 'GIT_CONFIG_NOSYSTEM', ...
        'OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', ...
        'OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_PATHDEF', ...
        'OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_STARTUP'};
    values = cellfun(@getenv, variables, 'UniformOutput', false);
    testCase.addTeardown(@() restore(state.root, old_path, old_dir, variables, values));
    setenv('HOME', state.home);
    setenv('XDG_CONFIG_HOME', fullfile(state.home, '.config'));
    setenv('GIT_CONFIG_NOSYSTEM', '1');
    setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_REGISTRY', state.registry);
    setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_PATHDEF', fullfile(state.root, 'pathdef.m'));
    setenv('OPTIPROFILER_MATLAB_PROBLEM_LIBRARY_STARTUP', fullfile(state.root, 'startup.m'));
    write_text(fullfile(state.home, '.gitconfig'), sprintf([ ...
        '[user]\n\tname = OptiProfiler Test\n\temail = test@example.org\n' ...
        '[url "file://%s/"]\n\tinsteadOf = https://github.com/\n' ...
        '[protocol "file"]\n\tallow = always\n' ...
        '[uploadpack]\n\tallowAnySHA1InWant = true\n' ...
        '[advice]\n\tdetachedHead = false\n'], state.remotes));

    test_dir = fileparts(mfilename('fullpath'));
    [~, attributes] = fileattrib(fullfile(test_dir, '..', '..', '..', '..'));
    repository = attributes.Name;
    src = fullfile(state.package, 'matlab', 'optiprofiler', 'src');
    mkdir(fileparts(src));
    copyfile(fullfile(repository, 'matlab', 'optiprofiler', 'src'), src);
    copyfile(fullfile(repository, 'setup.m'), state.package);
    copyfile(fullfile(repository, 'matlab_problem_libraries.lock'), state.package);
    s2 = fullfile(state.package, 'matlab', 'optiprofiler', 'problem_libs', 's2mpj');
    mkdir(fullfile(s2, 'src', 'matlab_problems'));
    write_text(fullfile(s2, 'README.md'), 'fixture');
end

function restore(root, old_path, old_dir, variables, values)
    path(old_path);
    cd(old_dir);
    clear setup
    for k = 1:numel(variables)
        setenv(variables{k}, values{k});
    end
    if isfolder(root)
        rmdir(root, 's');
    end
end

function output = run_setup(state, options)
    defaults = struct('install_matcutest', false, 'install_solar', false, ...
        'problem_library_root', state.libraries);
    for field = reshape(fieldnames(options), 1, [])
        defaults.(field{1}) = options.(field{1});
    end
    cd(state.package);
    clear setup
    output = evalc('setup(defaults)');
end

function files = solar_files()
    files = { ...
        'solar_select.m', sprintf('function names = solar_select(~)\nnames = {};\nend\n'); ...
        'solar_load.m', sprintf('function p = solar_load(~)\np = [];\nend\n'); ...
        'solar_collect_info.m', sprintf('function info = solar_collect_info()\ninfo = {};\nend\n')};
end

function path_name = remote_path(state, name)
    path_name = fullfile(state.remotes, 'optiprofiler', [name, '.git']);
end

function sha = make_remote(state, name, files)
    remote = remote_path(state, name);
    mkdir(remote);
    git(remote, 'init -q');
    git(remote, 'symbolic-ref HEAD refs/heads/main');
    for k = 1:size(files, 1)
        target = fullfile(remote, files{k, 1});
        if ~isfolder(fileparts(target)), mkdir(fileparts(target)); end
        write_text(target, files{k, 2});
    end
    git(remote, 'add -A');
    git(remote, 'commit -q -m base');
    sha = rev_parse(remote);
end

function sha = branch_commit(state, name)
% A commit that only a side branch reaches, as a pin on an unmerged branch is.
    remote = remote_path(state, name);
    git(remote, 'checkout -q -b pinned');
    write_text(fullfile(remote, 'pinned.txt'), 'pinned');
    git(remote, 'add -A');
    git(remote, 'commit -q -m pinned');
    sha = rev_parse(remote);
    git(remote, 'checkout -q main');
end

function clone(state, name, directory)
    target = fullfile(state.libraries, directory);
    git(state.root, sprintf('clone -q "https://github.com/optiprofiler/%s.git" "%s"', name, target));
end

function lock_commits(state, commits)
    lock_file = fullfile(state.package, 'matlab_problem_libraries.lock');
    text = fileread(lock_file);
    decoded = jsondecode(text);
    for name = reshape(fieldnames(commits), 1, [])
        for k = 1:numel(decoded.libraries)
            entry = decoded.libraries;
            if iscell(entry), entry = entry{k}; else, entry = entry(k); end
            if strcmp(entry.name, name{1})
                text = strrep(text, entry.commit, commits.(name{1}));
            end
        end
    end
    write_text(lock_file, text);
end

function sha = unknown_commit()
    sha = '0123456789abcdef0123456789abcdef01234567';
end

function sha = head(state, directory)
    sha = rev_parse(fullfile(state.libraries, directory));
end

function sha = rev_parse(directory)
    sha = strtrim(git(directory, 'rev-parse HEAD'));
end

function output = git(directory, arguments)
    [status, output] = system(sprintf('git -C "%s" %s', directory, arguments));
    assert(status == 0, 'git %s failed in %s: %s', arguments, directory, output);
end

function tf = is_registered(state, name)
    tf = false;
    if isfile(state.registry)
        data = load(state.registry);
        tf = any(strcmp({data.registrations.name}, name));
    end
end

function tf = on_path(directory)
    tf = any(strcmp(strsplit(path, pathsep), directory));
end

function write_text(filename, text)
    fid = fopen(filename, 'w');
    assert(fid >= 0, 'Cannot write %s.', filename);
    cleanup = onCleanup(@() fclose(fid));
    fwrite(fid, text, 'char');
end
