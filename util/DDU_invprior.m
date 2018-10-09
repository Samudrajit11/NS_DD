function [theta] = DDU_invprior(u,ranges)
%jeff=@(u,thmin,thmax) thmin*exp(u*log(thmax/thmin));
%theta(1)=tau
%theta(2)= sigma*sigma*tau
theta=zeros(length(u),1);
theta(1) = jeff(u(1),ranges(1,1),ranges(1,2));
theta(2) = jeff(u(2),ranges(2,1),ranges(2,2));


for i = 3:(length(u))
    
    theta(i) = u(i)*u(i);
    
end

end