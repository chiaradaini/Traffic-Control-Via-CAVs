% Work: Bus simulations
% Author: Maria Laura Delle Monache

% Velocity flow of cars


function[ve]=velocity(Rho,par)


%---------------------------------%
%      Data of the function       %
%---------------------------------%


%global Rho_max;
%global V_max;


%--------------------------------%
%      Body of the function      %
%--------------------------------%

ve=par.V_max*(1-Rho./par.Rho_max);
