%Plots a heatmap of wavefront intensity with desired settings.
function h = heatmap(x,y,wavefront)
    h = pcolor(x,y,wavefront); %plotting wavefront
    xlabel('x / m');
    ylabel('y / m');
    set(h,'edgecolor','none');
    colormap gray;
    axis square;
    colorbar('FontSize',24,'TickLabelInterpreter','latex');
end