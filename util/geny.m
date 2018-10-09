function [u]= geny(ranges,len)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generates an array u with the first and second 
%elements as uniformly distributed 
%random numbers and the rest of the elements
%as Ys which are Ornstein-Uhlenbeck process.
%len is defined as 'length of data + 1' in the skel file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
del_t=1;
u=zeros(len,1);
u(1)=rand;
u(2)=rand;
w=jeff(u(1),ranges(1,1),ranges(1,2));
z=jeff(u(2),ranges(2,1),ranges(2,2));

u(3) = sqrt(z) * erfinv((2*rand)-1);
p = exp(-del_t/w);
q = 1 - exp(-(2*del_t)/w);

for i = 4:(len)
    
    u(i) = (u(i-1)* p)+ (sqrt(z*q)* erfinv((2*rand)-1));
    
end

end





