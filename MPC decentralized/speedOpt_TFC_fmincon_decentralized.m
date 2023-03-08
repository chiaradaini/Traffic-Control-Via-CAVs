% Work: MPC implementation
% Author: Chiara Daini

%Optimisation of the speed based on the total fuel consumption

function [vel_opt] = speedOpt_TFC_fmincon_decentralized(Rho,par,v,vel,pos)

lb=30*ones(1,length(vel));      %lower limit
ub=100*ones(1,length(vel));     %upper limit

x0=vel;    %starting point

% options = optimset ('Algorithm','interior-point','TolX',1e-10,'TolFun',1e-10);
options = optimset ('Algorithm','interior-point','Display','iter','TolX',1e-10,'TolFun',1e-10);
% options = optimset ('Algorithm','sqp','Display','iter','TolX',1e-8,'TolFun',1e-5);

[vel_opt,~,~,~] = fmincon(@(vel)Optimizer_TFC_fmincon_decentralized(Rho,par,v,vel,pos),x0,[],[],[],[],lb,ub,[],options);


