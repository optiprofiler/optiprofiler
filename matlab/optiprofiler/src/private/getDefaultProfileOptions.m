function profile_options = getDefaultProfileOptions(solvers, feature, profile_options)

    if exist('parcluster', 'file') == 2 % Check if Parallel Computing Toolbox is available
        myCluster = parcluster();
        % Get the number of workers in the cluster
        nb_cores = myCluster.NumWorkers;
    else
        nb_cores = 1;
    end

    if ~isfield(profile_options, ProfileOptionKey.N_JOBS.value)
        profile_options.(ProfileOptionKey.N_JOBS.value) = nb_cores;
    end
    if ~isfield(profile_options, ProfileOptionKey.KEEP_POOL.value)
        profile_options.(ProfileOptionKey.KEEP_POOL.value) = true;
    end
    if ~isfield(profile_options, ProfileOptionKey.SEED.value)
        profile_options.(ProfileOptionKey.SEED.value) = 1;
    end
    if ~isfield(profile_options, ProfileOptionKey.BENCHMARK_ID.value)
        profile_options.(ProfileOptionKey.BENCHMARK_ID.value) = 'out';
    end
    if ~isfield(profile_options, ProfileOptionKey.SOLVER_NAMES.value)
        profile_options.(ProfileOptionKey.SOLVER_NAMES.value) = cellfun(@(s) func2str(s), solvers, 'UniformOutput', false);
        % Handle the case where the solver names are not valid MATLAB variable names.
        % Replace underscores with backslashes.
        profile_options.(ProfileOptionKey.SOLVER_NAMES.value) = cellfun(@(s) strrep(s, '_', '\_'), profile_options.(ProfileOptionKey.SOLVER_NAMES.value), 'UniformOutput', false);
    end
    if ~isfield(profile_options, ProfileOptionKey.SOLVER_ISRAND.value)
        profile_options.(ProfileOptionKey.SOLVER_ISRAND.value) = false(1, numel(solvers));
    end
    if ~isfield(profile_options, ProfileOptionKey.FEATURE_STAMP.value)
        profile_options.(ProfileOptionKey.FEATURE_STAMP.value) = getDefaultFeatureStamp(feature);
    end
    if ~isfield(profile_options, ProfileOptionKey.ERRORBAR_TYPE.value)
        profile_options.(ProfileOptionKey.ERRORBAR_TYPE.value) = 'minmax';
    end
    if ~isfield(profile_options, ProfileOptionKey.SAVEPATH.value)
        profile_options.(ProfileOptionKey.SAVEPATH.value) = pwd;
    end
    if ~isfield(profile_options, ProfileOptionKey.MAX_TOL_ORDER.value)
        profile_options.(ProfileOptionKey.MAX_TOL_ORDER.value) = 10;
    end
    if ~isfield(profile_options, ProfileOptionKey.MAX_EVAL_FACTOR.value)
        profile_options.(ProfileOptionKey.MAX_EVAL_FACTOR.value) = 500;
    end
    if ~isfield(profile_options, ProfileOptionKey.MERIT_FUN.value)
        profile_options.(ProfileOptionKey.MERIT_FUN.value) = @computeMeritValues;
    end
    if ~isfield(profile_options, ProfileOptionKey.PROJECT_X0.value)
        profile_options.(ProfileOptionKey.PROJECT_X0.value) = false;
    end
    if ~isfield(profile_options, ProfileOptionKey.RUN_PLAIN.value)
        profile_options.(ProfileOptionKey.RUN_PLAIN.value) = false;
    end
    if ~isfield(profile_options, ProfileOptionKey.SCORE_ONLY.value)
        profile_options.(ProfileOptionKey.SCORE_ONLY.value) = false;
    end
    if ~isfield(profile_options, ProfileOptionKey.SUMMARIZE_PERFORMANCE_PROFILES.value)
        profile_options.(ProfileOptionKey.SUMMARIZE_PERFORMANCE_PROFILES.value) = true;
    end
    if ~isfield(profile_options, ProfileOptionKey.SUMMARIZE_DATA_PROFILES.value)
        profile_options.(ProfileOptionKey.SUMMARIZE_DATA_PROFILES.value) = true;
    end
    if ~isfield(profile_options, ProfileOptionKey.SUMMARIZE_LOG_RATIO_PROFILES.value)
        profile_options.(ProfileOptionKey.SUMMARIZE_LOG_RATIO_PROFILES.value) = false;
    end
    if ~isfield(profile_options, ProfileOptionKey.SUMMARIZE_OUTPUT_BASED_PROFILES.value)
        profile_options.(ProfileOptionKey.SUMMARIZE_OUTPUT_BASED_PROFILES.value) = true;
    end
    if ~isfield(profile_options, ProfileOptionKey.SILENT.value)
        profile_options.(ProfileOptionKey.SILENT.value) = false;
    end
    if ~isfield(profile_options, ProfileOptionKey.SOLVER_VERBOSE.value)
        profile_options.(ProfileOptionKey.SOLVER_VERBOSE.value) = 1;
    end
    if ~isfield(profile_options, ProfileOptionKey.SEMILOGX.value)
        profile_options.(ProfileOptionKey.SEMILOGX.value) = true;
    end
    if ~isfield(profile_options, ProfileOptionKey.NORMALIZED_SCORES.value)
        profile_options.(ProfileOptionKey.NORMALIZED_SCORES.value) = false;
    end
    if ~isfield(profile_options, ProfileOptionKey.SCORE_WEIGHT_FUN.value)
        profile_options.(ProfileOptionKey.SCORE_WEIGHT_FUN.value) = @(x) 1;
    end
    if ~isfield(profile_options, ProfileOptionKey.SCORE_FUN.value)
        profile_options.(ProfileOptionKey.SCORE_FUN.value) = @(x) mean(x(:, :, 1, 1), 2);
    end
    if ~isfield(profile_options, ProfileOptionKey.LOAD.value)
        profile_options.(ProfileOptionKey.LOAD.value) = '';
    end
    if ~isfield(profile_options, ProfileOptionKey.LINE_COLORS.value)
        profile_options.(ProfileOptionKey.LINE_COLORS.value) = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.3010, 0.7450, 0.9330; 0.6350, 0.0780, 0.1840];
    end
    if ~isfield(profile_options, ProfileOptionKey.LINE_STYLES.value)
        profile_options.(ProfileOptionKey.LINE_STYLES.value) = {'-', '-.', ':', '--'};
    end
    if ~isfield(profile_options, ProfileOptionKey.LINE_WIDTHS.value)
        profile_options.(ProfileOptionKey.LINE_WIDTHS.value) = 1.5;
    end
    if ~isfield(profile_options, ProfileOptionKey.BAR_COLORS.value)
        % We will use the first two colors in the line color order as the bar colors.
        % If the line color order only has one color, we will use the same color for both bars.
        if iscell(profile_options.(ProfileOptionKey.LINE_COLORS.value))
            profile_options.(ProfileOptionKey.BAR_COLORS.value) = profile_options.(ProfileOptionKey.LINE_COLORS.value)(1:min(2, length(profile_options.(ProfileOptionKey.LINE_COLORS.value))));
        else
            profile_options.(ProfileOptionKey.BAR_COLORS.value) = profile_options.(ProfileOptionKey.LINE_COLORS.value)(1:min(2, size(profile_options.(ProfileOptionKey.LINE_COLORS.value), 1)), :);
        end
    end
end

function feature_stamp = getDefaultFeatureStamp(feature)
    % Generate a feature stamp to represent the feature with its options.

    switch feature.name
        case FeatureName.PERTURBED_X0.value
            % feature_name + perturbation_level + (distribution if it is gaussian or spherical)
            feature_stamp = sprintf('%s_%g', feature.name, feature.options.(FeatureOptionKey.PERTURBATION_LEVEL.value));
            if ischarstr(feature.options.(FeatureOptionKey.DISTRIBUTION.value)) && ismember(feature.options.(FeatureOptionKey.DISTRIBUTION.value), {'gaussian', 'spherical'})
                feature_stamp = sprintf('%s_%s', feature_stamp, feature.options.(FeatureOptionKey.DISTRIBUTION.value));
            end
        case FeatureName.NOISY.value
            % feature_name + noise_level + noise_type + (distribution if it is gaussian or uniform)
            feature_stamp = sprintf('%s_%g_%s', feature.name, feature.options.(FeatureOptionKey.NOISE_LEVEL.value), feature.options.(FeatureOptionKey.NOISE_TYPE.value));
            if ischarstr(feature.options.(FeatureOptionKey.DISTRIBUTION.value)) && ismember(feature.options.(FeatureOptionKey.DISTRIBUTION.value), {'gaussian', 'uniform'})
                feature_stamp = sprintf('%s_%s', feature_stamp, feature.options.(FeatureOptionKey.DISTRIBUTION.value));
            end
        case FeatureName.TRUNCATED.value
            % feature_name + significant_digits + (perturbed_trailing_zeros if it is true)
            feature_stamp = sprintf('%s_%d', feature.name, feature.options.(FeatureOptionKey.SIGNIFICANT_DIGITS.value));
            if feature.options.(FeatureOptionKey.PERTURBED_TRAILING_ZEROS.value)
                feature_stamp = sprintf('%s_perturbed_trailing_zeros', feature_stamp);
            end
        case FeatureName.LINEARLY_TRANSFORMED.value
            % feature_name + (rotated if it is true) + (condition_factor if it is not 0)
            feature_stamp = feature.name;
            if feature.options.(FeatureOptionKey.ROTATED.value)
                feature_stamp = sprintf('%s_rotated', feature_stamp);
            end
            if feature.options.(FeatureOptionKey.CONDITION_FACTOR.value) ~= 0
                feature_stamp = sprintf('%s_%g', feature_stamp, feature.options.(FeatureOptionKey.CONDITION_FACTOR.value));
            end
        case FeatureName.RANDOM_NAN.value
            % feature_name + nan_rate
            feature_stamp = sprintf('%s_%g', feature.name, feature.options.(FeatureOptionKey.NAN_RATE.value));
        case FeatureName.UNRELAXABLE_CONSTRAINTS.value
            % feature_name + (bounds if it is true) + (linear if it is true) + (nonlinear if it is true)
            feature_stamp = feature.name;
            if feature.options.(FeatureOptionKey.UNRELAXABLE_BOUNDS.value)
                feature_stamp = sprintf('%s_bounds', feature_stamp);
            end
            if feature.options.(FeatureOptionKey.UNRELAXABLE_LINEAR_CONSTRAINTS.value)
                feature_stamp = sprintf('%s_linear', feature_stamp);
            end
            if feature.options.(FeatureOptionKey.UNRELAXABLE_NONLINEAR_CONSTRAINTS.value)
                feature_stamp = sprintf('%s_nonlinear', feature_stamp);
            end
        case FeatureName.QUANTIZED.value
            % feature_name + mesh_size + (ground_truth if is_true it is true)
            feature_stamp = sprintf('%s_%g', feature.name, feature.options.(FeatureOptionKey.MESH_SIZE.value));
            if feature.options.(FeatureOptionKey.GROUND_TRUTH.value)
                feature_stamp = sprintf('%s_ground_truth', feature_stamp);
            end
        otherwise
            feature_stamp = feature.name;
    end
end