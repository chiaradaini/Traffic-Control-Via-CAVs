Rho(1,1:par.J/2)=0.3*par.Rho_max;                  
Rho(1,par.J/2+1:par.J)=0.3*par.Rho_max;

par.Rhol=Rho(1,1);      % For a simple Riemann Problem: \rho_L
par.Rhor=Rho(1,par.J);      % For a simple Riemann Problem: \rho_R