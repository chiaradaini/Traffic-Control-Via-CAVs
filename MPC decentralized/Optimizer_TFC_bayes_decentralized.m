% Work: MPC implementation
% Author: Chiara Daini

% Takes as input the vehicles and velocity to be optimized
% As output we have the total fuel consumption

function [TFC]=Optimizer_TFC_bayes_decentralized(Rho,par,v,vel,pos,xx)

    v{1}.vel=xx{1,1};
    v{1}.pos=pos(1);
    update(v{1}, par);

format long;
                     
dx(1,1:par.J)=par.L/par.J;                              % Grid size                   
                             
lambda=par.dt/par.dx_0;                             % Short-cut: Ratio time and space steps
x(1,1)=0;                                   % Starting point of the space grid
for i=2:par.J
    x(1,i)=x(1,i-1)+dx(1,i);                % Construction of the space grid
end

boundaries_check

    iter=0;
    T=par.hor_tot;
    n=length(v);
    m=[];
    t=zeros(1:iter);
    TFC=0;                                  %Total Fuel consumption
    ATT=0;                                  %Average Travel Time

    for i=1:n
        lane(i)=v{i}.lane;
        labels(i)=v{i}.label;
    end 
    l_tilde=unique(lane);

    for i=1:length(l_tilde)
            s(i,:)=(lane==l_tilde(i)*ones(size(lane)));
    end

    while T-par.hor_tot<=par.hor_opt 
        iter=iter+1;
        
        if iter==1
            t(iter)=par.hor_tot;
        else
            t(iter,1)=t(iter-1,1)+par.dt;
        end
        
        for i=1:n
                m(i)=v{i}.m;
        end 
        
        m_tilde=unique(m);
        
        clear z
        
        for i=1:length(m_tilde)
                z(i,:)=(m==m_tilde(i)*ones(size(m)));
        end
            
        for i=1:length(m_tilde)
            m_logic=find(z(i,:)==1); % trova l'indice dei veicoli nella stessa cella
            if length(m_logic) > 1   % hai n>1 veicoli nella stessa cella 
                clear s_tilde
                for j=1:length(m_logic)
                    s_tilde(:,j) = s(:,m_logic(j));
                end
                [height,width] = size(s_tilde);
                for p=1:height
                    s_logic=find(s_tilde(p,:)==1);
                    if length(s_logic)>1
                        clear A
                        clear B
                        for j=1:length(s_logic)
                            A(j,1)=v{m_logic(s_logic(j))}.vel;    % crea una matrice con le velocità
                            B(j,1)=v{m_logic(s_logic(j))}.pos;                        
                        end
                        for j=1:length(s_logic)
                            v{m_logic(s_logic(j))}.vel=min(A(A>0));
                            v{m_logic(s_logic(j))}.pos=max(B(B>0));
                            update(v{m_logic(s_logic(j))}, par);
                        end

                    else
                        %do nothing
                    end
                end
            else
                %do nothing
            end
        end
            
        % Godunov numerical flux
        num_flux = Gflux(Rho,par);
        %num_flux(1)=Fin;
        Dout = par.V_max*Rho(1,par.J).*(1-Rho(1,par.J)/par.Rho_max).*(Rho(1,par.J)<0.5*par.Rho_max) + 0.25*par.V_max*par.Rho_max.*(Rho(1,par.J)>=0.5*par.Rho_max);
        num_flux(par.J)=min(Fout,Dout);
    
    
        % Recontruction for classical shocks! When we have non-classical shocks they have precendence 
        RhoJ = Rho(1,1:par.J-2);
        RhoJ1 = Rho(1,3:par.J);
        index=(find((RhoJ < RhoJ1)));
        if ~isempty(index)
            for j=index(end):index(1)  % TO BE PARALLELIZED IF POSSIBLE
                k=j+1;
                if ~ismember(k,m)
                    l=(flux(Rho(1,k+1),par)-flux(Rho(1,k-1),par))/(Rho(1,k+1)-Rho(1,k-1));  %shock speed
                    dj=(Rho(1,k+1)-Rho(1,k))/(Rho(1,k+1)-Rho(1,k-1)); %reconstructed position in cell j
                    if l>0               
                        if dj>=0 && dj<=1
                            deltat=((1-dj)/l)*par.dx_0;
                            num_flux(k-1) = flux(Rho(1,k-1),par);
                            num_flux(k)=(min(deltat,par.dt)*flux(Rho(1,k+1),par)+max(par.dt-deltat,0)*flux(Rho(1,k-1),par))/par.dt;  
                        end
                    elseif l<0
                        if dj<=1 && dj>=0
                            deltatminus=((dj)/(-l))*par.dx_0; % what if l==0????
                            num_flux(k-1)=(min(deltatminus,par.dt)*flux(Rho(1,k-1),par)+max(par.dt-deltatminus,0)*flux(Rho(1,k+1),par))/par.dt;
                            num_flux(k) = flux(Rho(1,k+1),par);
                        end
                    end
                end
            end
        end
    
        % treatement of non-classical shocks TO BE PARALLELIZED IF POSSIBLE
        
        if (v{1}.m<par.J)
           if (v{1}.pos ~= x(1,v{1}.m)) 
                if (flux(RS(Rho(1,v{1}.m-1),Rho(1,v{1}.m+1),par,v{1}),par) > (v{1}.F_alpha+v{1}.vel*RS(Rho(1,v{1}.m-1),Rho(1,v{1}.m+1),par,v{1}))) 
                    v{1}.activity=1;
                    d = (v{1}.Rho_check - Rho(1,v{1}.m)) / (v{1}.Rho_check - v{1}.Rho_hat);
                    dtemp = par.dx_0*(1-d)/v{1}.vel ; 
                    if d>=0 && d<=1
                        num_flux(v{1}.m-1) = godunov(Rho(1,v{1}.m-1),v{1}.Rho_hat,par);
                        num_flux(v{1}.m) = (min(dtemp,par.dt)*flux(v{1}.Rho_check,par) + max(par.dt-dtemp,0)*flux(v{1}.Rho_hat,par))/par.dt;
                    end
                    v{1}.pos=v{1}.pos+v{1}.vel*par.dt;
                    v{1}.m=find(v{1}.pos>=x(1,:), 1, 'last' );
                elseif (flux(RS(Rho(1,v{1}.m-1),Rho(1,v{1}.m+1),par,v{1}),par) <= (v{1}.F_alpha+v{1}.vel*RS(Rho(1,v{1}.m-1),Rho(1,v{1}.m+1),par,v{1}))) && (Rho(1,v{1}.m-1)<Rho(1,v{1}.m+1))
                    v{1}.activity=0;
                    k=v{1}.m; % reconstruction of classical shock in cell m
                    l=(flux(Rho(1,k+1),par)-flux(Rho(1,k-1),par))/(Rho(1,k+1)-Rho(1,k-1));  %shock speed
                    dj=(Rho(1,k+1)-Rho(1,k))/(Rho(1,k+1)-Rho(1,k-1)); %reconstructed position in cell j
                    if l>0               
                        if dj>=0 && dj<=1
                            deltat=((1-dj)/l)*par.dx_0;
                            num_flux(k-1) = flux(Rho(1,k-1),par);
                            num_flux(k)=(min(deltat,par.dt)*flux(Rho(1,k+1),par)+max(par.dt-deltat,0)*flux(Rho(1,k-1),par))/par.dt;  
                        end
                    elseif l<0
                        if dj<=1 && dj>=0
                            deltatminus=((dj)/(-l))*par.dx_0; % what if l==0????
                            num_flux(k-1)=(min(deltatminus,par.dt)*flux(Rho(1,k-1),par)+max(par.dt-deltatminus,0)*flux(Rho(1,k+1),par))/par.dt;
                            num_flux(k) = flux(Rho(1,k+1),par);
                        end
                    end
                    v{1}.pos=v{1}.pos+Vbus(Rho(1,v{1}.m),par,v{1})*par.dt;
                    v{1}.m=find(v{1}.pos>=x(1,:), 1, 'last' );
                else  % nothing happens: no classical nor non-classical shocks
                    v{1}.activity=0;
                    v{1}.pos=v{1}.pos+Vbus(Rho(1,v{1}.m),par,v{1})*par.dt;
                    v{1}.m=find(v{1}.pos>=x(1,:), 1, 'last' );
                end
           else
                v{1}.activity=0;
                v{1}.pos=v{1}.pos+Vbus(Rho(1,v{1}.m),par,v{1})*par.dt;
                v{1}.m=find(v{1}.pos>=x(1,:), 1, 'last' );  
           end
        else                          
            v{1}.m=par.J;
            v{1}.pos=min(par.L,v{1}.pos+Vbus(Rho(1,v{1}.m),par,v{1})*par.dt) ;
        end

        yplot(1,iter)=v{1}.pos;

        % update of the flux
        Fin = min(fin(T,par), par.V_max*Rho(1,1).*(1-Rho(1,1)/par.Rho_max).*(Rho(1,1)>0.5*par.Rho_max) + 0.25*par.V_max*par.Rho_max.*(Rho(1,1)<=0.5*par.Rho_max) );
        num_fluxm = [Fin num_flux(1,1:end-1)];
        Rho(1,:)=Rho(1,:)-par.lambda*(num_flux(1,:)-num_fluxm(1,:)); %density update
        
        %Total Fuel Consumption
        fc=fuel_consumption(Rho,par);
        TFC = TFC + sum(fc) * par.dx_0 * par.dt;
        
        %Average travel time
        tt=AverageTravelTime(Rho,par);
        ATT=ATT+sum(tt) * par.dx_0 * par.dt;

        Rhoplot(iter,1:par.J)=Rho;
        T=T+par.dt;
        
    end
    
end
