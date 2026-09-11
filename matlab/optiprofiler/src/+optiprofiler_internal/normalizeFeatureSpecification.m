function specification = normalizeFeatureSpecification(input, varargin)
%NORMALIZEFEATURESPECIFICATION One pure, callback-free normalization boundary.
% Preserve inherited numeric 0/1 logicals, selected NaN/Inf acceptance, and the
% absence of perturbation-level validation; a policy change is separate work.
    flat = parseFlatOptions(varargin);
    if isText(input)
        tokens = strsplit(char(input),'+','CollapseDelimiters',false);
        entries = cell(1,numel(tokens)); route = 'feature_name'; accepted = {};
        for k = 1:numel(tokens)
            name = atomicName(tokens{k});
            definition = optiprofiler_internal.featureDefinitions(name);
            routed = struct(); fields = fieldnames(flat);
            for j = 1:numel(fields)
                if ismember(fields{j},definition.local_keys)
                    routed.(fields{j}) = flat.(fields{j}); accepted{end+1} = fields{j};
                end
            end
            entries{k} = struct('name',name,'options',routed);
        end
        if ~all(ismember(fieldnames(flat),accepted))
            error('MATLAB:Feature:InvalidOptionForFeature', 'A supplied option is not accepted by any declared stage.');
        end
    elseif isa(input,'struct') && isscalar(input)
        entries = {input}; route = 'feature';
    elseif isa(input,'cell') && ~isempty(input)
        entries = reshape(input,1,[]); route = 'feature';
    else
        error('MATLAB:Feature:FeaturenameNotString', ...
            'Feature requires a name, scalar stage struct, nonempty cell array or canonical Feature.');
    end
    if strcmp(route,'feature') && ~isempty(fieldnames(flat))
        error('MATLAB:Feature:StructuredOverrides', 'Structured input keeps local options inside each stage, not as extra arguments.');
    end
    stages = {}; declared_entries = {}; occurrences = struct(); names = {}; stochastic = false;
    for k = 1:numel(entries)
        entry = entries{k}; supplied = struct();
        if isa(entry,'struct') && isscalar(entry) && isfield(entry,'name')
            name = entry.name;
            if isfield(entry,'options'), supplied = entry.options; end
            if ~all(ismember(fieldnames(entry),{'name','options'}))
                error('MATLAB:Feature:InvalidStage', 'A stage accepts only name and options.');
            end
        elseif isText(entry)
            name = entry;
        else
            error('MATLAB:Feature:InvalidStage', 'Each stage is scalar text or a scalar struct containing name and optional options.');
        end
        name = atomicName(name);
        definition = optiprofiler_internal.featureDefinitions(name);
        supplied = optionStruct(supplied); options = definition.defaults; fields = fieldnames(supplied);
        for j = 1:numel(fields)
            key = fields{j};
            if ~ismember(key,definition.local_keys)
                error('MATLAB:Feature:InvalidOptionForFeature', 'Option %s is not local to stage %s.',key,name);
            end
            options.(key) = validateOption(name,key,supplied.(key));
        end
        declared_entries{end+1} = struct('name',name,'options',supplied); names{end+1} = name;
        if strcmp(name,'plain'), continue; end
        if ~isfield(occurrences,name), occurrences.(name) = 0; end
        occurrence = occurrences.(name); occurrences.(name) = occurrence+1;
        stages{end+1} = struct('name',name,'occurrence',occurrence, ...
            'identity',sprintf('%s#%d',name,occurrence),'code',definition.code,'options',options);
        definition = optiprofiler_internal.featureDefinitions(name,options);
        stochastic = stochastic || definition.is_stochastic;
    end
    effective_name = 'plain';
    if ~isempty(stages), effective_name = strjoin(cellfun(@(s) s.name,stages,'UniformOutput',false),'+'); end
    declaration = struct('route',route,'entries',{declared_entries});
    specification = struct('stages',{stages},'declared',declaration,'name',effective_name, ...
        'declared_name',strjoin(names,'+'),'is_stochastic',stochastic,'is_identity',isempty(stages), ...
        'specification_version','matlab-feature-spec-v2');
end

function flat = parseFlatOptions(values)
    if isempty(values), flat = struct(); return; end
    if numel(values)==1 && isa(values{1},'struct'), flat = optionStruct(values{1}); return; end
    if mod(numel(values),2)~=0
        error('MATLAB:Feature:InvalidNumberOfArguments', 'Feature options must be a scalar struct or name/value pairs.');
    end
    flat = struct();
    for k = 1:2:numel(values)
        key = values{k};
        if ~isText(key) || ~isvarname(char(key))
            error('MATLAB:Feature:UnknownOption', 'Option names must be scalar text naming a valid field.');
        end
        flat.(lower(char(key))) = values{k+1};
    end
    flat = optionStruct(flat);
end

function output = optionStruct(input)
    if ~isa(input,'struct') || ~isscalar(input)
        error('MATLAB:Feature:InvalidStage', 'Stage options must be a scalar struct.');
    end
    output = struct(); fields = fieldnames(input);
    definitions = optiprofiler_internal.featureDefinitions(); known = [definitions.local_keys];
    for k = 1:numel(fields)
        key = lower(fields{k});
        if strcmp(key,'n_runs')
            error('MATLAB:Feature:ExperimentOption', ...
                'n_runs is an experiment option. Use benchmark(solvers, struct(''n_runs'', N, ...)), not Feature(..., ''n_runs'', N).');
        end
        if ~ismember(key,known), error('MATLAB:Feature:UnknownOption', 'Unknown option for feature: %s.',key); end
        output.(key) = input.(fields{k});
    end
end

function name = atomicName(value)
    if ~isText(value), error('MATLAB:Feature:InvalidStage', 'Stage names must be scalar text.'); end
    name = lower(strtrim(char(value)));
    if isempty(name) || contains(name,'+')
        error('MATLAB:Feature:InvalidStage', 'Stage names must be nonempty and atomic; use one entry per stage.');
    end
end

function value = validateOption(name,key,value)
    switch key
        case 'distribution'
            if ~isa(value,'function_handle') && ~isText(value), fail('distribution_NotFunctionHandle','distribution must be text or a function handle.'); end
            if isText(value)
                value = char(value);
                if strcmp(name,'noisy') && ~ismember(value,{'gaussian','uniform'})
                    fail('distribution_NotFunctionHandle_noisy','noisy distribution must be gaussian, uniform or a function handle.');
                elseif strcmp(name,'perturbed_x0') && ~ismember(value,{'gaussian','spherical'})
                    fail('distribution_NotFunctionHandle_perturbed_x0','perturbed_x0 distribution must be gaussian, spherical or a function handle.');
                end
            end
        case 'nan_rate'
            if ~realScalar(value) || value<0 || value>1, fail('nan_rate_NotBetween_0_1','nan_rate must be between zero and one.'); end
        case 'significant_digits'
            if ~realScalar(value) || rem(value,1)~=0 || value<=0, fail('significant_digits_NotPositiveInteger','significant_digits must be a positive integer.'); end
        case 'noise_level'
            if ~realScalar(value) || value<0, fail('noise_level_NotPositive','noise_level must be nonnegative.'); end
        case 'condition_factor'
            if ~(realScalar(value) && value>=0), fail('condition_factor_InvalidInput','condition_factor must be nonnegative.'); end
        case 'mesh_size'
            if ~realScalar(value) || value<=0, fail('mesh_size_NotPositive','mesh_size must be positive.'); end
        case 'noise_type'
            if ~isText(value) || ~ismember(char(value),{'absolute','relative','mixed'}), fail('noise_type_InvalidInput','noise_type must be absolute, relative or mixed.'); end
            value = char(value);
        case 'noise_mode'
            if ~isText(value) || ~ismember(char(value),{'random','deterministic'}), fail('noise_mode_InvalidInput','noise_mode must be random or deterministic.'); end
            value = char(value);
        case 'mesh_type'
            if ~isText(value) || ~ismember(char(value),{'absolute','relative'}), fail('mesh_type_InvalidInput','mesh_type must be absolute or relative.'); end
            value = char(value);
        case 'noise_map'
            if ~isa(value,'function_handle') && ~isText(value), fail('noise_map_NotFunctionHandle','noise_map must be text or a function handle.'); end
            if isText(value)
                value = char(value);
                if ~strcmp(value,'chebyshev'), fail('noise_map_InvalidInput','The named noise_map must be chebyshev.'); end
            end
        case {'perturbed_trailing_digits','rotated','unrelaxable_bounds','unrelaxable_linear_constraints','unrelaxable_nonlinear_constraints','ground_truth'}
            if ~logicalScalar(value), fail([key,'_NotLogical'],'The option must be a logical scalar or numeric zero/one.'); end
        case {'mod_x0','mod_bounds','mod_linear_ub','mod_linear_eq','mod_affine','mod_fun','mod_cub','mod_ceq'}
            if ~isa(value,'function_handle'), fail([key,'_NotFunctionHandle'],'The modifier must be a function handle.'); end
        case 'perturbation_level'
            % Preserve the historical absence of an extra value validator.
    end
end

function value = isText(input)
    value = (isa(input,'char') && (isrow(input) || isempty(input))) || ...
        (isa(input,'string') && isscalar(input) && ~ismissing(input));
end
function value = realScalar(input)
    value = builtin('isnumeric',input) && isreal(input) && isscalar(input);
end
function value = logicalScalar(input)
    value = (isa(input,'logical') && isscalar(input)) || (realScalar(input) && (input==0 || input==1));
end
function fail(suffix,message)
    error(['MATLAB:Feature:',suffix],'%s',message);
end
