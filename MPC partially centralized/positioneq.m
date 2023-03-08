% Work: Bus simulations
% Author: Maria Laura Delle Monache

% Flux function
% Concave flux function

function[pp]=positioneq(t,x,lambda,time,bart,barx,par)


%---------------------------------%
%      Data of the function       %
%---------------------------------%





%--------------------------------%
%      Body of the function      %
%--------------------------------%

pp=(position(t,x,time,bart,barx,par)-x)/t-lambda;