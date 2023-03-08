% Work: Bus simulations
% Author:

% Inflow of cars


function [Fin]=fin(t,par)


%---------------------------------%
%      Data of the function       %
%---------------------------------%


%global f_max;
%global Tf;
%global Rhol;


%--------------------------------%
%      Body of the function      %
%--------------------------------%
% 
    if t<0.5*par.Tf
         Fin = par.f_max;
    else
        Fin = 0.;
    end
    
