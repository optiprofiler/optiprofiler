function [record, reason] = propagateProblemReference(record, name, options, reason)
%PROPAGATEPROBLEMREFERENCE Retain raw target facts or the first fallback reason.
%   Built-ins retain scoring truth, including original-coordinate residuals
%   under affine changes. Optional x_ref remains source metadata in original
%   coordinates, with no transport, oracle query, or random-stream draw.
%   Ground-truth quantization, every custom stage, and unknown stages use the
%   empirical target. An already invalidated record stays invalidated, and
%   its first reason is carried through later stages unchanged.

    if nargin < 4, reason = ''; end
    if isempty(record)
        record = [];
        return
    end
    switch char(name)
        case {'plain', 'noisy', 'truncated', 'random_nan', 'nonquantifiable_constraints', ...
                'unrelaxable_constraints', 'perturbed_x0', 'permuted', 'linearly_transformed'}
            retained = true;
        case 'quantized'
            % Retained only for an explicit logical false; a missing or
            % non-logical value counts as ground truth.
            retained = isstruct(options) && isscalar(options) && isfield(options, 'ground_truth') ...
                && islogical(options.ground_truth) && isscalar(options.ground_truth) && ~options.ground_truth;
        case 'custom'
            % This release deliberately has no custom reference callback.
            retained = false;
        otherwise
            retained = false;
    end
    % Scoring residuals remain in original coordinates even if the solver's
    % constraints use a transformed representation. Keep nonzero maxcv too.
    % x_ref is original-coordinate metadata, so it needs no matrix operation.
    if ~retained
        record = [];
        if strcmp(char(name),'custom'), reason = 'custom_feature';
        elseif strcmp(char(name),'quantized'), reason = 'quantized_ground_truth';
        else, reason = 'unsupported_feature'; end
    else
        reason = '';
    end
end
