% Work: MPC with multiple CAVs with respect of Total fuel consumption
% Author: Chiara Daini, Paola Goatin

% MPC centralized: the vehicles velocity is optimized on the base of the
% knoledge of the density along the whole highway.
% The optimisation is made with fmincon and it's made at the same moment
% for each vehicle

% Optimizaton made only with bayesOpt + fmincon

clc
clear all

tic

%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%

par.V_max = 140;
par.Rho_max = 400;
par.L = 50;
par.J = par.L / 0.2;
par.f_max = 0.25 * par.V_max * par.Rho_max;
par.dx_0 = par.L /par.J;
par.cfl = 0.9 / par.V_max;
par.dt = par.cfl * par.dx_0;
par.lambda=par.dt/par.dx_0;

% MPC one step horizon
par.Tf=1;   %time horizon = 1 hour
par.hor_opt=65*par.dt;   %time horizon for optimisation = 5 min
par.hor_avant=52*par.dt;    %time horizon for simulation before restarting optimization
par.hor_apres=13*par.dt;    %remaining time horizon after restarted optimization
par.hor=65*par.dt; %time horizon for simulation = 5 min
par.hor_tot=0;  %time horizon for the while cycle

% MPC receding horizon
% par.Tf=1;   %time horizon = 1 hour
% par.hor_opt=0.25;   %time horizon for optimisation = 15 min
% par.hor=65*par.dt; %time horizon for simulation = 5 min
% par.hor_tot=0;  %time horizon for the while cycle

% 1 hour optimization
% par.Tf=1; %time horizon = 1 hour
% par.hor_opt=1; %time horizon for optimisation = 1 min
% par.hor=1; %time horizon for simulation = 5 min
% par.hor_tot=0; %time horizon for the while cycle

par.t_0=0;
TFC=0;
yplot_opt=[];
Rhoplot_opt=[];
velplot=[];
t=[];

initial_datum

%%%%%%%%%%%%%%%% VEICHLES %%%%%%%%%%%%%%%%%

vel(1)=50;
pos(1)=4.5;
v{1}=av(pos(1),1,0.6,vel(1),1); 
update(v{1}, par);

vel(2)=50;
pos(2)=9;
v{2}=av(pos(2),2,0.6,vel(2),2);
update(v{2}, par);

% vel(3)=50;
% pos(3)=13.5;
% v{3}=av(pos(3),3,0.6,vel(3),3);
% update(v{3}, par);
% 
% vel(4)=50;
% pos(4)=18;
% v{4}=av(pos(4),1,0.6,vel(4),4);
% update(v{4}, par);
% 
% vel(5)=50;
% pos(5)=22.5;
% v{5}=av(pos(5),2,0.6,vel(5),5);
% update(v{5}, par);
% 
% vel(6)=50;
% pos(6)=27;
% v{6}=av(pos(6),3,0.6,vel(6),6);
% update(v{6}, par);
% 
% vel(7)=50;
% pos(7)=31.5;
% v{7}=av(pos(7),1,0.6,vel(7),7);
% update(v{7}, par);
% 
% vel(8)=50;
% pos(8)=36;
% v{8}=av(pos(8),2,0.6,vel(8),8);
% update(v{8}, par);
% 
% vel(9)=50;
% pos(9)=41.5;
% v{9}=av(pos(9),3,0.6,vel(9),9);
% update(v{9}, par);
% 
% vel(10)=50;
% pos(10)=46;
% v{10}=av(pos(10),1,0.6,vel(10),10);
% update(v{10}, par);

%%%%%%%%%%%%%%%%% MPC %%%%%%%%%%%%%%%%%%%

while par.hor_tot < par.Tf
    
    if par.hor_tot == 0 %first optimization
        [vel_opt] = speedOpt_TFC_bayes_centralized(Rho, par, v, vel, pos);

        vel_bayes = vel_opt.XAtMinObjective{1,:};

        [vel_opt_fmincon] = speedOpt_TFC_fmincon_centralized(Rho, par, v, vel_bayes, pos);

        vel = vel_opt_fmincon;

        if par.hor_tot + par.hor > par.Tf
            par.hor = par.Tf - par.hor_tot;
        end
        
        par.hor = par.hor_avant; %set par.hor equal the first period of time
        
        [cost_opt, pos_new, yplot, Rhoplot, Rho_new, t_new]=Simulator_TFC(Rho, par, v, vel, pos); %simulation for the first period of time

        Rho = Rho_new; %update the density

        pos = pos_new; %update the positions of the vehicles
        
        t=[t; t_new]; %update the time
        
        TFC= TFC + cost_opt; %update the TFC
        
        par.hor_tot = par.hor_tot + par.hor;

        yplot_opt = [yplot_opt yplot];
        Rhoplot_opt = [Rhoplot_opt; Rhoplot];
        velplot = [velplot; vel];
        
        %the following optimization starts before the end of the simulation
        
        [vel_opt] = speedOpt_TFC_bayes_centralized(Rho, par, v, vel, pos); 

        vel_bayes = vel_opt.XAtMinObjective{1,:};

        [vel_opt_fmincon] = speedOpt_TFC_fmincon_centralized(Rho, par, v, vel_bayes, pos);

        vel2 = vel_opt_fmincon;
        
        par.hor = par.hor_apres; %set par.hor equal the second period of time
        
        [cost_opt, pos_new, yplot, Rhoplot, Rho_new, t_new]=Simulator_TFC(Rho, par, v, vel, pos); %simulation for the second period of time
                
        par.hor_tot = par.hor_tot + par.hor;
        
        Rho = Rho_new; %update the density

        pos = pos_new; %update the positions of the vehicles
        
        t=[t; t_new]; %update the time
        
        TFC= TFC + cost_opt; %update the TFC

        yplot_opt = [yplot_opt yplot];
        Rhoplot_opt = [Rhoplot_opt; Rhoplot];
        velplot = [velplot; vel];
        
        vel = vel2;
        
    else %all the other optimization

        if par.hor_tot + par.hor > par.Tf
            par.hor = par.Tf - par.hor_tot;
        end
        
        par.hor = par.hor_avant; %set par.hor equal the first period of time
        
        [cost_opt, pos_new, yplot, Rhoplot, Rho_new, t_new]=Simulator_TFC(Rho, par, v, vel, pos); %simulation for the first period of time

        Rho = Rho_new; %update the density

        pos = pos_new; %update the positions of the vehicles
        
        t=[t; t_new]; %update the time
        
        TFC= TFC + cost_opt; %update the TFC
        
        par.hor_tot = par.hor_tot + par.hor;
        
        yplot_opt = [yplot_opt yplot];
        Rhoplot_opt = [Rhoplot_opt; Rhoplot];
        velplot = [velplot; vel];
        
        %the following optimization starts before the end of the simulation        
        
        [vel_opt] = speedOpt_TFC_bayes_centralized(Rho, par, v, vel, pos);

        vel_bayes = vel_opt.XAtMinObjective{1,:};

        [vel_opt_fmincon] = speedOpt_TFC_fmincon_centralized(Rho, par, v, vel_bayes, pos);

        vel2 = vel_opt_fmincon;
        
        par.hor = par.hor_apres; %set par.hor equal the second period of time
        
        [cost_opt, pos_new, yplot, Rhoplot, Rho_new, t_new]=Simulator_TFC(Rho, par, v, vel, pos); %simulation for the second period of time
                
        par.hor_tot = par.hor_tot + par.hor;
        
        Rho = Rho_new; %update the density

        pos = pos_new; %update the positions of the vehicles
        
        t=[t; t_new]; %update the time
        
        TFC= TFC + cost_opt; %update the TFC

        yplot_opt = [yplot_opt yplot];
        Rhoplot_opt = [Rhoplot_opt; Rhoplot];
        velplot = [velplot; vel];
        
        vel = vel2;
   
    end

end 

%%%%%%%%%%%%%%%%%%%%% SAVE DATA %%%%%%%%%%%%%%%%%%%%%%
% 
% i=length(v);
% 
% density = sprintf('TFC_MPC_Rhoplotcentral_duplex_%d', i);
% vehicleposition = sprintf('TFC_MPC_yplotcentral_duplex_%d', i);
% speed = sprintf('TFC_MPC_velplotcental_duplex_%d', i);
% 
% save(density,'Rhoplot_opt')
% save(vehicleposition,'yplot_opt')
% save(speed,'velplot')

toc




