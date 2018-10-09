function logl = DD_logl_m(obs,theta)
%computes log-likelihood

T = length(obs(1,:));  
D = theta(3:(T+2));
del_t=1;
logl = 0; 
dim=2; %dimension
for t=1:T
    if (D(t)==0)
        D(t)=10^(-30);
    end    
    logl = logl - (0.5*dim*log(4*pi*D(t)*del_t)) - (((obs(1,t)*obs(1,t))+(obs(2,t)*obs(2,t)))/(4 * D(t)*del_t));
end    


  
end



