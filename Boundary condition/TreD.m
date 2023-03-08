load('Rhoplot')
load('yplot')

% function[]=TreD(Rhoplot2)

% time=linspace(0,5,max(size(RhoMatrix)));
% space=1:size(RhoMatrix,1);
figure;
map = multigradient([0 0.8 0; 1 1 0; 0.8 0 0]); 
colormap(gca,map);
% imagesc(time,1:7,RhoMatrix);
imagesc(x1,t,Rhoplot);
hold on
plot(yplot,t,'k','linewidth',2)
set(gca,'YDir','normal');

shading flat

ylabel 'time [h]';
xlabel 'space [h]';
title 'Density';

% title 'Penetration rate';
%  export_fig(['DensityH','.pdf'],'-pdf','-transparent','-p8');

hold off

% saveas(gca,'different_lanes_figure4_a.png');
