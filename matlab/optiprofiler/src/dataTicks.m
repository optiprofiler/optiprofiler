function [ticks, tickLabels] = dataTicks(ratio_cut_data)

    if ratio_cut_data >= 5
        max_power = floor(ratio_cut_data);
        ticks = linspace(0, max_power, 6);
        ticks(2:end-1) = round(ticks(2:end-1));
        ticks = unique(ticks, 'stable');
    elseif ratio_cut_data >= 1
        max_power = floor(ratio_cut_data);
        ticks = (0:1:max_power);
    elseif ratio_cut_data >= 1e-1
        ticks = [0 ratio_cut_data];
    else
        ticks = [0];
    end

    tickLabels = arrayfun(@(x) num2str(2 ^ x - 1), ticks, 'UniformOutput', false);
end