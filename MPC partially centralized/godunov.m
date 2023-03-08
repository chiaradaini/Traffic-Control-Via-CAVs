% Work: bus simulations
% Author: Maria Laura Delle Monache

% Godunov scheme

function[god]=godunov(u,v,par)

%---------------------------------------%
%         Data of the function          %
%---------------------------------------%

%global Rho_max

%---------------------------------------%
%         Body of the function          %
%---------------------------------------%


if u<=v
    god=min(flux(u,par),flux(v,par));
elseif (v<u) && (u<=0.5*par.Rho_max)
    god=flux(u,par);
elseif (v<0.5*par.Rho_max) && (0.5<u*par.Rho_max)
    god=flux(0.5*par.Rho_max,par);
else
    god=flux(v,par);
end

