function test()
    
    fig = figure('visible', 'on');
    ax = axes(fig);
    % Create a 3-D tensor that follows linear decay in the third dimension.
    y = zeros(8, 5, 100);
    for i = 1:100
        y(:, :, i) = -i;
    end
    % Multiply with a random matrix.
    y = y .* rand(8, 5, 100);


    labels = {'a', 'b', 'c','d', 'e', 'f','g','h'};
    ismaxcv = false;
    profile_options = struct('range_type', 'minmax');
    drawHist(ax, y, labels, ismaxcv, profile_options)

    
end