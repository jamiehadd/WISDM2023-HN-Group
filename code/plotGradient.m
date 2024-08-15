function fig = plotGradient(X1,Y1,X2,Y2,step_size)
    clc; close all
    fig = figure;
    semilogy(median(X1),  median(Y1), [ '-', 'r'],'LineWidth',2)
    hold on
    semilogy(median(X2),  median(Y2), [ '-.', 'b'],'LineWidth',2)
    
    for conf = 0:step_size:(1-step_size)
        upper_bound = quantile(Y1, conf + step_size);
        x_upper_bound = quantile(X1, conf + step_size);
        lower_bound = quantile(Y1,conf);
        x_lower_bound = quantile(X1,conf);
        fill([x_upper_bound fliplr(x_lower_bound)], [upper_bound fliplr(lower_bound)], [conf 0 0], 'FaceAlpha', 0.2, 'LineStyle', 'none')
    end
    
    for conf = 0:step_size:(1-step_size)
        upper_bound = quantile(Y2, conf + step_size);
        x_upper_bound = quantile(X2, conf + step_size);
        lower_bound = quantile(Y2,conf);
        x_lower_bound = quantile(X2,conf);
        fill([x_upper_bound fliplr(x_lower_bound)], [upper_bound fliplr(lower_bound)], [0 0 conf], 'FaceAlpha', 0.2, 'LineStyle', 'none')
    end
    
    %colorbar work
    color_interval = 0:step_size:(1-step_size);
    map = zeros(length(color_interval),3); 
    for row = 1:length(color_interval)
        old_rgb = [color_interval(row), color_interval(row), color_interval(row)];
        new_rgb = 1 - 0.2.*(1 - old_rgb);
        map(row,:) = new_rgb;
    end

    colormap(map);
    c = colorbar;
    c.Label.String = "quantile";
    c.Label.FontSize = 14;
    c.Label.Interpreter = 'latex';
end