%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Bus v2.0 - ML DELLE MONACHE       %
%            Bus velocity subroutine         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[vb]=Vbus(Rho,par,av)

%global V_b;

vb = min(av.vel,vel(Rho,par));
    