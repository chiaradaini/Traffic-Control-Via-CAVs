function [ FC_rate ] =fuel_consumption(Rho,par)
%Fuel consumption function
%   Function that calculates the consumed fuel depending on the velocity

s=par.V_max*(1-Rho./par.Rho_max);
K=5.7*10^(-12)*(s.^6) -3.6*10^(-9)*(s.^5)+7.6*10^(-7)*(s.^4) -6.1*10^(-5)*(s.^3) +1.9*10^(-3)*(s.^2) +1.6*10^(-2)*s+0.99;    % s is in km/h, k liters/h

FC_rate = K.*Rho;

end


