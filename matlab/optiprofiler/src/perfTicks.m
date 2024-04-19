function [ticks, tickLabels] = perfTicks(ratio_cut_perf)

    if ratio_cut_perf >= 5
        max_power = floor(ratio_cut_perf);
        ticks = linspace(0, max_power, 6);
        ticks(2:end-1) = round(ticks(2:end-1));
        ticks = unique(ticks, 'stable');
    elseif ratio_cut_perf >= 1
        max_power = floor(ratio_cut_perf);
        ticks = (0:1:max_power);
    elseif ratio_cut_perf >= 1e-3
        ticks = [0 ratio_cut_perf];
    else
        ticks = [0];
    end

    tickLabels = arrayfun(@(x) num2str(2 ^ x), ticks, 'UniformOutput', false);
end