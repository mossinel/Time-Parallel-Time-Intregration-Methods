function [t,u_delta]=Euler_Delta_theta(df,f,tspan,u0,theta,check)
dt=(tspan(2)-tspan(1));
% since it is only related to the coarse integrator, it is sufficient to
% define only one coarse step of integration
[t,u]=theta_nonlinear(df,f,tspan, u0, 1,theta);
u_delta=eye(length(u0),length(u0));
Dphi_B=@(tt,u) eye(length(u0))-dt*df(tt,u)*(1-theta);
Dphi_F=@(tt,u) eye(length(u0))+dt*df(tt,u)*(theta);
F=Dphi_F(t(1),u(1,:));
G=Dphi_B(t(end),u(end,:));
if strcmp(check,"Y")         
                func=@(tt,uu) uu+dt*theta*ff(tt,uu);
                [flag, difference]=check_jacobian(F,func,1,u(1,:),tspan(1),tspan(end));
                if max(difference)>=1e-5
                    error("!!wrong explicit jacobian!!")
                end
                [flag, difference]=check_jacobian_BE(G,f,df,1,u,tspan(1),tspan(2),theta);
                if max(difference)>=1e-5
                    error("!!wrong implicit jacobian!!")
                end
end
            
u_delta=G\(F*u_delta); 
  
end
