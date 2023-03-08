% THIS FILE CAN PRINT CURVES FOR VELOCITIES OR DENSITY STARTING FORM THE
% DATA INCLUDED IN THE FOLDER 'DATA'

par.L = 50;
par.J = par.L / 0.2;
par.Rho_max = 400;
[m,n]=size(Rhoplot_opt);
t = linspace(0,1,m);

dx(1,1:par.J)=par.L/par.J;                              % Grid size 
x(1,1)=0;                                   % Starting point of the space grid
for i=2:par.J
    x(1,i)=x(1,i-1)+dx(1,i);                % Construction of the space grid
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SPEEDS PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stairs(velplot, 'linewidth', 1.5);
title('Optimal speed');
xlabel('Time [h]','fontsize',15);
ylabel('speed [km/h]','fontsize',15);
ylim([20, 110]);
xlim([0,25]);
%%%legend('CAV 1','CAV 2','CAV 3','CAV 4','CAV 5','CAV 6','CAV 7','CAV 8','CAV 9','CAV 10');

%%%%%%%%%%%%%%%%%%%%%%%%%%% DENSITY PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% imagesc(x,t,Rhoplot_opt);
% map = multigradient([0 0.8 0; 1 1 0; 0.8 0 0]);
% colormap(map);
% c=colorbar;
% c.Limits=[0,par.Rho_max];
% hold on
% plot(yplot_opt,t,'k','linewidth',2)
% set(gca,'YDir','normal');
% 
% shading flat
% 
% ylabel('time [h]','fontsize',20);
% xlabel('space [km]', 'fontsize',20);
% %title('Density','fontsize',20);
% 
% hold off

