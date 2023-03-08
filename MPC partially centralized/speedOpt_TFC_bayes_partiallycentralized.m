% Work: MPC implementation
% Author: Chiara Daini

%Optimisation of the speed based on the total fuel consumption

function [vel_opt] = speedOpt_TFC_bayes_partiallycentralized(Rho,par,v,vel,pos, v_neigh, vel_neigh, pos_neigh)

% concatenate x + 1 for the first vehicle
tmp = char("x"+1);

% upper and lower limits for the first vehicle
lim = optimizableVariable(tmp,[30,100],'Type','real');

% upper and lower limits for the other vehicles
if length(v) > 1
    for i = 2:length(v)  
        lim = [lim optimizableVariable(char("x"+i),[30,100],'Type','real')];
    end
end

xx = vel;

fun = @(xx)Optimizer_TFC_bayes_partiallycentralized(Rho,par,v,vel,pos,xx, v_neigh, vel_neigh, pos_neigh);

max_eval = 15 * length(v);

% vel_opt = bayesopt(fun,lim,'AcquisitionFunctionName','expected-improvement-per-second-plus','ExplorationRatio',2);
%  vel_opt = bayesopt(fun,lim,'PlotFcn', [],'Verbose',0,'MaxObjectiveEvaluations',20);
vel_opt = bayesopt(fun,lim,'PlotFcn', [], 'MaxObjectiveEvaluations',max_eval, 'UseParallel', true);
% vel_opt = bayesopt(fun,lim,'AcquisitionFunctionName','probability-of-improvement');

% PlotFcn', [] = not plotting the curves
% 'Verbose',0 = not displaying the results
% 'MaxObjectiveEvaluations',20 = max number of iterations
