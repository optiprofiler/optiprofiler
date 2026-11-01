function receipt = prepareMatcutestRuntime(adapter_root, runtime_parent)
%PREPAREMATCUTESTRUNTIME validates an adapter's exact runtime path receipt.
%   The supported platforms are those of the matcutest entry of the lock.
    platforms = lockedPlatforms();
    if ~ismember(matlabPlatformName(), platforms)
        error('OptiProfiler:MatcutestUnsupportedPlatform', ...
            'MatCUTEst setup supports only the platforms of its entry in matlab_problem_libraries.lock (%s), not %s.', ...
            strjoin(platforms, ', '), matlabPlatformName());
    end
    entrypoint = fullfile(adapter_root, 'matcutest_setup.m');
    if ~isfile(entrypoint)
        error('OptiProfiler:MatcutestReceiptRequired', ...
            'Adapter %s has no managed runtime entry point. Adopt the reviewed adapter fix; native setup is not run.', adapter_root);
    end
    original_path = path;
    original_dir = pwd;
    cleanup = onCleanup(@() restore_context(original_path, original_dir)); %#ok<NASGU>
    cd(adapter_root);
    addpath(adapter_root);
    clear matcutest_setup
    receipt = matcutest_setup(runtime_parent);
    fields = {'schema_version', 'provider', 'runtime_root', 'runtime_paths'};
    if ~isstruct(receipt) || ~isscalar(receipt) || ~all(isfield(receipt, fields)) || ...
            ~isequal(receipt.schema_version, 1) || ~isequal(receipt.provider, 'matcutest') || ...
            ~ischar(receipt.runtime_root) || ~iscellstr(receipt.runtime_paths) || ...
            ~isrow(receipt.runtime_paths) || isempty(receipt.runtime_paths)
        error('OptiProfiler:InvalidMatcutestReceipt', 'Invalid MatCUTEst runtime receipt.');
    end
    expected_root = fullfile(canonical_directory(runtime_parent), 'matcutest');
    expected_path = fullfile(expected_root, 'mtools', 'src');
    if ~strcmp(canonical_directory(expected_root), expected_root) || ...
            ~strcmp(receipt.runtime_root, expected_root) || ...
            ~isequal(receipt.runtime_paths, {expected_path}) || ...
            ~strcmp(canonical_directory(expected_path), expected_path)
        error('OptiProfiler:InvalidMatcutestReceipt', ...
            'MatCUTEst receipt must name only its canonical managed mtools/src directory.');
    end
end

function platforms = lockedPlatforms()
    % A lock that cannot say which platforms are supported says nothing about
    % this platform.
    try
        readMatlabProblemLibraryLock();
    catch cause
        error('OptiProfiler:MatcutestLockUnavailable', ...
            'MatCUTEst setup reads its platforms from matlab_problem_libraries.lock, which cannot be read: %s', cause.message);
    end
    entry = lockedProblemLibraryEntry('matcutest');
    if ~isstruct(entry) || ~isfield(entry, 'platforms') || isempty(entry.platforms)
        error('OptiProfiler:MatcutestLockUnavailable', ...
            'matlab_problem_libraries.lock lists no platforms for MatCUTEst, so MatCUTEst setup cannot check this platform.');
    end
    platforms = reshape(cellstr(entry.platforms), 1, []);
end

function directory = canonical_directory(directory)
    quoted = ['''', strrep(char(directory), '''', '''"''"'''), ''''];
    [status, output] = system(['realpath -e -- ', quoted]);
    if status ~= 0 || ~isfolder(strtrim(output))
        error('OptiProfiler:InvalidMatcutestReceipt', 'Invalid runtime directory: %s', directory);
    end
    directory = strtrim(output);
end

function restore_context(old_path, old_dir)
    path(old_path);
    cd(old_dir);
end
