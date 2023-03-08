% Work: MPC with multiple CAVs with respect of Total fuel consumption
% Author: Chiara Daini, Paola Goatin

% MPC decentralized: the vehicles velocity is optimized on the base of the
% knoledge of the density along the whole highway.
% The optimisation is made with fmincon and it's made separately for each
% vehicle

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

OPT = 0.25;
SIM = 65*par.dt;

% MPC
par.Tf=1;   %time horizon = 1 hour
par.hor_opt=OPT;   %time horizon for optimisation = 15 min
par.hor=SIM; %time horizon for simulation = 5 min
par.hor_tot=0;  %time horizon for the while cycle

% 1 hour optimization
% par.Tf=1; %time horizon = 1 hour
% par.hor_opt=1; %time horizon for optimisation = 1 min
% par.hor=1; %time horizon for simulation = 5 min
% par.hor_tot=0; %time horizon for the while cycle

par.ray = 5;

%%%%%%%%%%%%%%%% VEICHLES %%%%%%%%%%%%%%%%%

V0 = 50;

vel(1)=V0;
pos(1)=4.5;
v{1}=av(pos(1),1,0.6,vel(1),1); 
update(v{1}, par);

vel(2)=V0;
pos(2)=9;
v{2}=av(pos(2),2,0.6,vel(2),2);
update(v{2}, par);

vel(3)=V0;
pos(3)=13.5;
v{3}=av(pos(3),3,0.6,vel(3),3);
update(v{3}, par);

vel(4)=V0;
pos(4)=18;
v{4}=av(pos(4),1,0.6,vel(4),4);
update(v{4}, par);

vel(5)=V0;
pos(5)=22.5;
v{5}=av(pos(5),2,0.6,vel(5),5);
update(v{5}, par);

vel(6)=V0;
pos(6)=27;
v{6}=av(pos(6),3,0.6,vel(6),6);
update(v{6}, par);

vel(7)=V0;
pos(7)=31.5;
v{7}=av(pos(7),1,0.6,vel(7),7);
update(v{7}, par);

vel(8)=V0;
pos(8)=36;
v{8}=av(pos(8),2,0.6,vel(8),8);
update(v{8}, par);

vel(9)=V0;
pos(9)=41.5;
v{9}=av(pos(9),3,0.6,vel(9),9);
update(v{9}, par);

vel(10)=V0;
pos(10)=46;
v{10}=av(pos(10),1,0.6,vel(10),10);
update(v{10}, par);


%%%%%%%%%%%%%%%%% MPC %%%%%%%%%%%%%%%%%%%

for i = 1:length(v)
    
    par.t_0=0;
    TFC=0;
    yplot_opt=[];
    Rhoplot_opt=[];
    velplot=[];
    t=[];
    initial_datum
    
    v_loop=v(1,1:i);
    vel_loop=vel(1,1:i);
    pos_loop=pos(1,1:i);

    vel_init = vel_loop;
    
    for i = 1:length(vel_init)
        labels(i) = v{i}.label;
    end
    
    while par.hor_tot < par.Tf
        
        pos_step = pos_loop;
        vel_step = vel_init;
        %cycle for the decentralisation: each vehicle is analyzed by itslef
        %considering the density of the street but not the presence of the
        %other CAVs
        for j = 1:length(vel_loop)
            
            v_array = v_loop(1,j);
            vel_array = vel_init(1,j);
            pos_array = pos_loop(1,j);
            
            lpos = v_array{1}.pos - par.ray;
            upos = v_array{1}.pos + par.ray;   
            
            neigh = 0;
            neigh_labels = 0;
    
            for p = 1:length(v_loop)
                if p == j
                    % do nothing
                elseif v_loop{p}.pos >= lpos && v_loop{p}.pos <= upos
                        neigh_labels(p) = v_loop{p}.label;
                else
                    % do nothing
                end         
            end
            neigh_labels = neigh_labels(neigh_labels>0);
            
            for p = 1:length(neigh_labels)
                neigh(p) = find(neigh_labels(p) == labels);
            end
            
            if neigh == 0
    
            % Optimisation through fmincon for 15 min
    
                vel_opt = speedOpt_TFC_bayes_decentralized(Rho, par, v_array, vel_array, pos_array);
    
                vel_bayes = vel_opt.XAtMinObjective{1,:};
    
                vel_opt_fmincon = speedOpt_TFC_fmincon_decentralized(Rho, par, v_array, vel_bayes, pos_array);
    
                vel_init(j) = vel_opt_fmincon;
            
            else
                neigh = neigh(neigh>0); 
    
                for p = 1:length(neigh)
                    vel_neigh(p) = v{neigh(p)}.vel;
                    pos_neigh(p) = v{neigh(p)}.pos;
                    v_neigh(p) = v(1,neigh(p));
                end
    
                vel_group = [vel_array vel_neigh];
                pos_group = [pos_array pos_neigh];
                v_group = [v_array v_neigh];

                % Optimisation through fmincon for 15 min

                [vel_opt] = speedOpt_TFC_bayes_centralized(Rho, par, v_group, vel_group, pos_group);

                vel_bayes = vel_opt.XAtMinObjective{1,:};

                vel_opt_fmincon = speedOpt_TFC_fmincon_centralized(Rho, par, v_group, vel_bayes, pos_group);

                vel_fmincon = vel_opt_fmincon;

                vel_init(j) = vel_fmincon(1);

                clear vel_opt
                clear vel_bayes
                clear vel_opt_fmincon
                clear vel_fmincon

    
            end
            
             clear vel_neigh
             clear pos_neigh
             clear v_neigh
             
             clear vel_group
             clear pos_group
             clear v_group 
             
            for k = 1:length(vel_init)
                v_loop{k}.vel = vel_step(k);
                v_loop{k}.pos = pos_step(k);
                update(v{k},par);
            end
    
        end
        
        if par.hor_tot + par.hor > par.Tf
            par.hor = par.Tf - par.hor_tot;
        end
        
        % Simulation is made with all the velocity obtained from the previous
    %     optimisation
        % Simulator for 5 min
        [cost_opt, pos_new, yplot, Rhoplot, Rho_new, t_new]=Simulator_TFC(Rho, par, v_loop, vel_init, pos_loop);
        
        Rho = Rho_new;
        
        pos_loop = pos_new;
    
        par.hor_tot = par.hor_tot + par.hor;
        
        t=[t; t_new];
        
        TFC = TFC + cost_opt;
            
        yplot_opt = [yplot_opt yplot];
        Rhoplot_opt = [Rhoplot_opt; Rhoplot];
        velplot = [velplot; vel_init];
    
    end 


    TFC_array(i) = TFC;
    %reset the parameters
    if par.hor_opt == OPT
        par.hor=SIM;  %MPC  
    else
        par.hor=1;   % 1 hour
    end
    par.hor_tot=0;

        %%%%%%%%%%%%%%%%%%%%% SAVE DATA %%%%%%%%%%%%%%%%%%%%%%

    % MPC
    density = sprintf('MPC_Rhoplot_quasidecentral_%dCAV', i);
    vehicleposition = sprintf('MPC_yplot_quasidecentral_%dCAV', i);
    speed = sprintf('MPC_velplot_quasidecentral_%dCAV', i);
    
%     % 1 hour
%     density = sprintf('TFC_1hour_Rhoplotdecentral_duplex_%d', i);
%     vehicleposition = sprintf('TFC_1hour_yplotdecentral_duplex_%d', i);
%     speed = sprintf('TFC_1hour_velplotdecental_duplex_%d', i);

    save(density,'Rhoplot_opt')
    save(vehicleposition,'yplot_opt')
    save(speed,'velplot')
% 
end
% 
save('MPC_quasidecentral_TFC_array', 'TFC_array');

toc




