% Work: Bus simulations
% Author: Maria Laura Delle Monache

% Flux function
% Concave flux function

function[f]=flux(Rho,par)


%---------------------------------%
%      Data of the function       %
%---------------------------------%


%global Rho_max;
%global V_max;


%--------------------------------%
%      Body of the function      %
%--------------------------------%

if (Rho>0) && (Rho<par.Rho_max)
    f=par.V_max*Rho*(1-Rho/par.Rho_max);
else
    f=0;
end