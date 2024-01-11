function testplot
    
    styles = {'-','-o', '-^'};
    colors = [1 0 0; 0 0 1; 0 1 0];
    ax = axes(); 
    ax.LineStyleOrder = styles; 
    ax.ColorOrder = colors;
    ax.LineStyleCyclingMethod = 'withcolor'; 
    hold on
    plot(0:.1:1, rand(1,11)+(1:6)')
    legend
end