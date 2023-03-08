classdef av < handle
    
    properties
        pos; %position of the vehicle (to be assigned)
        lane; %lane of the vehicle (to be assigned)
        lane_occupied; %number of laned occupied by one vehicle. normally is one, by two vehicles are coupled could grow up to 3 (to be assigned)
        vel; %speed of the vehicle (to be assigned)
        label; %label of the vehicle (to be assigned)
        alpha; %portion of street occupied by the vehicle %calculated automatically
        rhoalpha; %calculated automatically
        f_alpha; %calculated automatically
        F_alpha; %calculated automatically
        Rho_check; %calculated automatically
        Rho_hat; %calculated automatically
        Rho_star; %calculated automatically
        Rho_cr; %calculated automatically
        m; % cell occupied by the vehicle (calculated automatically)
        activity; % activity of the vehicle (1 = active; 0 = not active)(calculated automatically)
    end
    
     
    methods
        function obj = av(pos, lane, lane_occupied, vel, label)
            if nargin == 5
                obj.pos = pos;
                obj.lane = lane;
                obj.lane_occupied = lane_occupied;
                obj.vel = vel;
                obj.label = label;
                obj.alpha = NaN;                
                obj.rhoalpha = NaN;
                obj.f_alpha = NaN;
                obj.F_alpha = NaN;
                obj.Rho_check = NaN;
                obj.Rho_hat=NaN;
                obj.Rho_star=NaN;
                obj.Rho_cr=NaN;
                obj.m=NaN;
                obj.activity=0;
            end
        end

        function obj = updatealpha(obj, par)
            obj.alpha = (par.lane_tot - obj.lane_occupied) / par.lane_tot;
        end

        function  out=myavproperties(obj, par)
            out.rhoalpha = 0.5*[obj.alpha]*[par.Rho_max]*(1-[obj.vel]/[par.V_max]);
            if obj.alpha==0.
                out.f_alpha = 0;
            else
                out.f_alpha = [par.V_max]*[out.rhoalpha]*(1-[out.rhoalpha]/([obj.alpha]*[par.Rho_max]));
            end 
            out.F_alpha = [out.f_alpha] - [obj.vel]*[out.rhoalpha];
           
            out.Rho_check = (0.5/[par.V_max])*([par.Rho_max]*([par.V_max]-[obj.vel]) - sqrt(([par.Rho_max]*([par.V_max]-[obj.vel]))^2 -4*[par.V_max]*[par.Rho_max]*[out.F_alpha]));
            
            out.Rho_hat = (0.5/[par.V_max])*([par.Rho_max]*([par.V_max]-[obj.vel]) + sqrt(([par.Rho_max]*([par.V_max]-[obj.vel]))^2 -4*[par.V_max]*[par.Rho_max]*[out.F_alpha]));
 
            out.Rho_star = [par.Rho_max]*(1-[obj.vel]/[par.V_max]);
           
            out.Rho_cr = 0.5*[par.Rho_max]*(1-[obj.vel]/[par.V_max]);
            
            out.m = round(obj.pos/par.dx_0)+1;
           
        end
        
        function obj = update(obj, par)
            out=myavproperties(obj, par);
            obj.rhoalpha = out.rhoalpha;
            obj.f_alpha = out.f_alpha;
            obj.F_alpha =out.F_alpha;
            obj.Rho_check =out.Rho_check;
            obj.Rho_hat = out.Rho_hat;
            obj.Rho_star = out.Rho_star;
            obj.Rho_cr = out.Rho_cr;
            obj.m = out.m;
        end
    end
            
end