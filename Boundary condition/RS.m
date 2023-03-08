% Work:Bus 3.0
% Author: Maria Laura Delle Monache

% Classical Riemann Solver at x/t=V_b

function[r]=RS(Rhol,Rhor,par,av)


%---------------------------------%
%      Data of the function       %
%---------------------------------%


%global Rho_cr;
%global V_b;
%global V_max;
%global Rho_max;


%--------------------------------%
%      Body of the function      %
%--------------------------------%

lambda=(flux(Rhor,par)-flux(Rhol,par))/(Rhor-Rhol);
V=av.vel;
    
% Shock
if Rhol<=Rhor
    
    if V >= lambda
        r=Rhor;
    else
        r=Rhol;
    end
% Rarefaction waves    
else
    if V >= par.V_max*(1-2*Rhor/par.Rho_max)
        r=Rhor;
    elseif V <= par.V_max*(1-2*Rhol/par.Rho_max)
        r=Rhol;
    else
        r=av.Rho_cr;
    end
end