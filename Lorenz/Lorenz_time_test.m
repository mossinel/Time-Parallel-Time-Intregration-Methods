function []= Lorenz_time_test(sigma,r,b,T,MF,MG,N,u0,K)


for t=T
    [U,u,err_inf,err_2]=Lorentz_Parareal_new(sigma, r , b,t,MF,MG,N,u0,K);
    figure
    semilogy([0:K],err_inf,'Linewidth',2)
    hold on
    
    semilogy([0:K],err_2,'Linewidth',2)
    legend("norm L^{\infty}","L^2",'Fontsize',15)
    title(["error, T="+T])
    xlabel('iteration')
    set(gca, 'FontSize', 15);
    
    figure
    subplot(3,1,1)
    plot([0:T/(length(u)-1):T],u(:,1),'r',[0:T/(length(u)-1):T],u(:,2),'g',[0:T/(length(u)-1):T],u(:,3), ...
    'b',[0:T/(length(U{1})-1):T],U{1}(:,1),'ro-',[0:T/(length(U{1})-1):T],U{1}(:,2),'go-',[0:T/(length(U{1})-1):T],U{1}(:,3),'bo-','Linewidth',2)
    xlabel('time') 
    ylabel('x')
    set(gca, 'FontSize', 10);

    subplot(3,1,2)
    plot([0:T/(length(u)-1):T],u(:,1),'r',[0:T/(length(u)-1):T],u(:,2),'g',[0:T/(length(u)-1):T],u(:,3), ...
        'b',[0:T/(length(U{K/2})-1):T],U{K/2}(:,1),'ro-',[0:T/(length(U{K/2})-1):T],U{K/2}(:,2),'go-',[0:T/(length(U{K/2})-1):T],U{K/2}(:,3),'bo-','Linewidth',2)
     xlabel('time') 
    ylabel('x')
    set(gca, 'FontSize', 10);

    subplot(3,1,3)
    plot([0:T/(length(u)-1):T],u(:,1),'r',[0:T/(length(u)-1):T],u(:,2),'g',[0:T/(length(u)-1):T],u(:,3), ...
    'b',[0:T/(length(U{K})-1):T],U{K}(:,1),'ro-',[0:T/(length(U{K})-1):T],U{K}(:,2),'go-',[0:T/(length(U{K})-1):T],U{K}(:,3),'bo-','Linewidth',2)
     xlabel('time') 
    ylabel('x')
    set(gca, 'FontSize', 10);
end