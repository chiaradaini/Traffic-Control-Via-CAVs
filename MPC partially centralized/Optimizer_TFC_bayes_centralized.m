% Work: MPC implementation
% Author: Chiara Daini

% Takes as input the vehicles and velocity to be optimized
% As output we have the total fuel consumption

function [TFC]=Optimizer_TFC_bayes_centralized(Rho,par,v,vel,pos,xx)

for i=1:length(v)
    v{i}.vel=xx{1,i};
    v{i}.pos=pos(i);
    update(v{i}, par);
end

format long;
                     
dx(1,1:par.J)=par.L/par.J;                              % Grid size                   
                             
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

        % In this lines, the matter is to verify if two vehicle in the
        % same cell have the same line too. In this case the vehicle
        % behind will follow the first one and adapt its own velocity
        % and position

        for i=1:length(m_tilde)
            m_logic=find(z(i,:)==1); % trova l'indice dei veicoli nella stessa cella
            if length(m_logic) > 1   % hai n>1 veicoli nella stessa cella 
                clear s_tilde
                for j=1:length(m_logic)
                    s_tilde(:,j) = s(:,m_logic(j));
                end
                [height,~] = size(s_tilde);
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
%                     else
%                         %do nothing
                    end
                end
%             else
%                 %do nothing
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
      
        % In the following lines, the matter is to verify which vehicle is
        % the active and to create a list of vehicle to be treated in the
        % for loop. This is important to avoid oscillation of the density
        % during the simulation
        active1=[];
        notactive1=[];
        active2 = [];
        active3 = [];
        notactive2 = [];
        notactive3 = [];
        active=[];
        notactive=[];
        pos_act=[];
        pos_notact=[];
        
        for i=1:length(v)
            if v{i}.activity==1
                active1(i) = v{i}.label;
            elseif v{i}.activity==0
                notactive1(i) = v{i}.label;
            end
        end 
        
        active2 = active1(active1>0);
        for k = 1:length(active2)
            active3(k) = find(active2(k) == labels);
            pos_act(k) = v{active3(k)}.pos;
        end
        [~,idx_act] = sort(pos_act,'descend');
        active = active3(idx_act);

        notactive2 = notactive1(notactive1>0);
        for k = 1:length(notactive2)
            notactive3(k) = find(notactive2(k) == labels);
            pos_notact(k) = v{notactive3(k)}.pos;
        end
        [~,idx_act] = sort(pos_notact,'descend');
        notactive = notactive3(idx_act);
        
        act=[notactive active];
          
        for p = 1:n

        if (v{act(p)}.m<par.J)
           if (v{act(p)}.pos ~= x(1,v{act(p)}.m)) 
                if (flux(RS(Rho(1,v{act(p)}.m-1),Rho(1,v{act(p)}.m+1),par,v{act(p)}),par) > (v{act(p)}.F_alpha+v{act(p)}.vel*RS(Rho(1,v{act(p)}.m-1),Rho(1,v{act(p)}.m+1),par,v{act(p)}))) 
                    v{act(p)}.activity=1;
                    d = (v{act(p)}.Rho_check - Rho(1,v{act(p)}.m)) / (v{act(p)}.Rho_check - v{act(p)}.Rho_hat);
                    dtemp = par.dx_0*(1-d)/v{act(p)}.vel ; 
                    if d>=0 && d<=1
                        num_flux(v{act(p)}.m-1) = godunov(Rho(1,v{act(p)}.m-1),v{act(p)}.Rho_hat,par);
                        num_flux(v{act(p)}.m) = (min(dtemp,par.dt)*flux(v{act(p)}.Rho_check,par) + max(par.dt-dtemp,0)*flux(v{act(p)}.Rho_hat,par))/par.dt;
                    end
                    v{act(p)}.pos=v{act(p)}.pos+v{act(p)}.vel*par.dt;
                    v{act(p)}.m=find(v{act(p)}.pos>=x(1,:), 1, 'last' );
                elseif (flux(RS(Rho(1,v{act(p)}.m-1),Rho(1,v{act(p)}.m+1),par,v{act(p)}),par) <= (v{act(p)}.F_alpha+v{act(p)}.vel*RS(Rho(1,v{act(p)}.m-1),Rho(1,v{act(p)}.m+1),par,v{act(p)}))) && (Rho(1,v{act(p)}.m-1)<Rho(1,v{act(p)}.m+1))
                    v{act(p)}.activity=0;
                    k=v{act(p)}.m; % reconstruction of classical shock in cell m
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
                    v{act(p)}.pos=v{act(p)}.pos+Vbus(Rho(1,v{act(p)}.m),par,v{act(p)})*par.dt;
                    v{act(p)}.m=find(v{act(p)}.pos>=x(1,:), 1, 'last' );
                else  % nothing happens: no classical nor non-classical shocks
                    v{act(p)}.activity=0;
                    v{act(p)}.pos=v{act(p)}.pos+Vbus(Rho(1,v{act(p)}.m),par,v{act(p)})*par.dt;
                    v{act(p)}.m=find(v{act(p)}.pos>=x(1,:), 1, 'last' );
                end
           else
                v{act(p)}.activity=0;
                v{act(p)}.pos=v{act(p)}.pos+Vbus(Rho(1,v{act(p)}.m),par,v{act(p)})*par.dt;
                v{act(p)}.m=find(v{act(p)}.pos>=x(1,:), 1, 'last' );  
           end
        else                          
            v{act(p)}.m=par.J;
            v{act(p)}.pos=min(par.L,v{act(p)}.pos+Vbus(Rho(1,v{act(p)}.m),par,v{act(p)})*par.dt) ;
        end

        yplot(act(p),iter)=v{act(p)}.pos;

        end
    
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
