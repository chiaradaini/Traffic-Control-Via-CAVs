% Work: Bus simulations
% Author: Maria Laura Delle Monache

% Flux function
% Concave flux function

function[p]=position(t,x,time,bart,barx,par)


%---------------------------------%
%      Data of the function       %
%---------------------------------%

%global V_max;



%--------------------------------%
%      Body of the function      %
%--------------------------------%

p=x+par.V_max*(t-time)+sqrt(t-time)*((barx-x-par.V_max*bart)/sqrt(bart));