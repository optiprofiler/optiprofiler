function y_shift = computeHistoryYShift(history, profile_options)
%COMPUTEHISTORYYSHIFT Shared translation for native and portable history plots.
%   HISTORY must already be a finite display copy from processHistYaxes.

    y_shift = 0;
    if strcmp(profile_options.(ProfileOptionKey.ERRORBAR_TYPE.value), 'meanstd')
        y_mean = squeeze(mean(history, 2));
        % Match prepareHistoryPlotData's established MATLAB sample bands.
        y_std = squeeze(std(history, 0, 2));
        y_lower = y_mean - y_std;
        y_min = min(y_lower(:));
    else
        y_min = min(history(:));
    end

    % Positive tiny values are valid on log axes; adding eps would erase
    % their relative variation. Shift only genuinely nonpositive data.
    if any(diff(history(:))) && y_min <= 0
        y_shift = max(eps - y_min, eps(-y_min) - y_min);
    end
end
