function F = Gflux(Rho,par)
%global V_max;
%global Rho_max;
% demand
D = par.V_max*Rho.*(1-Rho/par.Rho_max).*(Rho<0.5*par.Rho_max) + 0.25*par.V_max*par.Rho_max.*(Rho>=0.5*par.Rho_max);
% supply
S = par.V_max*Rho.*(1-Rho/par.Rho_max).*(Rho>0.5*par.Rho_max) + 0.25*par.V_max*par.Rho_max.*(Rho<=0.5*par.Rho_max);
S = [S(2:end) S(end)];
% Godunov flux
F = min(D,S);
end  