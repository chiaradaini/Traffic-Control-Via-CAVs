function [ inv_v ] = AverageTravelTime(Rho,par)
%Average Travel Time

%Function that returns a vector containing the inverse of the velocity at every time sample   

s=par.V_max*(1-Rho./par.Rho_max);

inv_v = 1./ s;






