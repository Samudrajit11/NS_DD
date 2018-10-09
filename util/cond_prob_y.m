function [prob] = cond_prob_y(u)
%calculates conditional probability of {Y} given \tau and D_*. 
% len is the number of Ys. 
%   \tau is u(1), D_* is u(2) and the rest of the elements of u are Ys.
del_t=1;
len = length(u)-2;
z = zeros((len-1),1);
p = 1-exp(-(2*del_t)/u(1)) ;
r = exp(-del_t/u(1)) ;
x = ((pi*u(2))^(-len/2))* (p)^(-(len-1)/2);
y = u(2)* p;
for k=1:(len-1)
    z(k)=(u(k+3)-u(k+2)* r)^(2);
end    
S = sum(z); 
prob = x * exp(-S/y)*exp(-u(3)*u(3)/u(2));


end

