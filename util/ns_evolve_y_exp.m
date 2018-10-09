function [walker_new,step_mod,stp_mod]=ns_evolve_y_exp(obs,model,logLstar,walker,step_mod,stp_mod,ranges)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uses the MCMC technique to find a new sample uniformly distributed inside
% the region of parameter space with higher likelihood than the minimum 
% requirement (logLstar). 
%
% Some of the arguments of the function are
% 
% obs - a 2xT matrix of observations
% walker - the walker that constitutes the starting point of the
%   MCMC process.
% logLstar - the minimum requirement for the likelihood of the new sample
% step_mod - a scalar value that regulates the average step length of the
%   subsequent call of the function for tau and D_*. When the remaning parameter space
%   becomes small, the MCMC steps are adjusted in length to ensure a
%   success rate of the MCMC steps of about 50%.
%stp_mod - step modifier for the Ys.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize step_mode if run for the first time
if step_mod==0
  step_mod=1;   
end

if stp_mod==0
    stp_mod=1;
end    
Sm = stp_mod ;
lent = length(walker.u);

%Specify likelihood function
logl = model.logl;

% Attempt to generate new walker through independent guess
walker_new.u = model.genu();
walker_new.theta = model.invprior(walker_new.u);
walker_new.logl=logl(obs,walker_new.theta); % Calculates new likelihood 

if(walker_new.logl <= logLstar)	% Do MCMC if likelihood requirement failed
   % Find a new walker via a MCMC process.
   % Return input values and logLstar if no steps succeed
   walker_new = walker;
   walker_new.logl=logLstar;
  


   i = 0;
   reject1 = 0; %counts failed steps for tau and D_* 
   accept2 = 0; %counts succesful steps for Ys
   reject2 = 0; %counts failed steps for Ys
   old=zeros(lent,1);
   new=zeros(lent,1);
   
   while(i < model.options.nsteps)
      % Propose step for parameters
      % Ensure that new walker.u is between 0 and 1
      
      
     for k=1:2 
         walker.u(k) = mod(walker_new.u(k) +(rand-0.5)* step_mod,1);
        
     end
      walker.theta(1) = jeff(walker.u(1),ranges(1,1),ranges(1,2));
      walker.theta(2) = jeff(walker.u(2),ranges(2,1),ranges(2,2));
      new(1:2)=walker.theta(1:2);
      new(3:lent)=walker.u(3:lent);
      old(1:2)=walker_new.theta(1:2);
      old(3:lent)=walker_new.u(3:lent);
      prob_old=cond_prob_y(old);
      prob_new=cond_prob_y(new);
      prob_rat = prob_new/prob_old ;
      
      
      check = rand;
         
      %walker.logl=logl(obs,walker.theta);
      
      
      %if(walker.logl<=logLstar)
          
         
          %for k=1:2
              %walker.u(k)=walker_new.u(k);
              %walker.theta(k)=walker_new.theta(k);
          %end
          %reject1 = reject1 + 1 ;
      if(prob_rat < 1 && check>prob_rat)
           for k=1:2
              walker.u(k)=walker_new.u(k);
              walker.theta(k)=walker_new.theta(k);
          end
          reject1 = reject1 + 1 ;
      end
       
      
      
      
      Gamma1= walker.u(4)*exp(-1/walker.theta(1));
      sigsq1 = 0.5 * walker.theta(2)*(1-exp(-2/walker.theta(1)));
      
      
      walker.u(3)= Gamma1 + (( walker.u(3)-Gamma1)*sqrt(1-(Sm*Sm))) + (Sm*sqrt(sigsq1)*randn);
      walker.theta(3) = walker.u(3) * walker.u(3) ;
      
      walker.logl=logl(obs,walker.theta);
      
      if(walker.logl > logLstar)
         
          accept2 = accept2 + 1 ;
      else
          walker.u(3)=walker_new.u(3);
          walker.theta(3)=walker_new.theta(3);
          reject2 = reject2 + 1 ;
      end
      walker.logl=logl(obs,walker.theta);
      sigsq = 0.5 * walker.theta(2) * tanh(1/walker.theta(1));      
      for j=4:(lent-1)
          
            Gamma =0.5*(walker.u(j-1)+walker.u(j+1))*(1/cosh(1/walker.theta(1)));
          
            walker.u(j)= Gamma + ((walker.u(j)-Gamma)*sqrt(1-(Sm*Sm))) + (Sm*sqrt(sigsq)*randn);
            walker.theta(j)= walker.u(j) * walker.u(j);
            
            F=obs(j-1)-obs(j-2);
            
            walker.logl= walker.logl-log(sqrt(4*pi*walker.theta(j))) - ((F*F)/(4 * walker.theta(j)))+log(sqrt(4*pi*walker_new.theta(j))) + ((F*F)/(4 * walker_new.theta(j))); 
            
            if(walker.logl > logLstar)
                
                accept2 = accept2 + 1 ;
            else
		walker.logl= walker.logl+log(sqrt(4*pi*walker.theta(j))) +((F*F)/(4 * walker.theta(j)))-log(sqrt(4*pi*walker_new.theta(j))) -((F*F)/(4 * walker_new.theta(j))); 
                walker.u(j)=walker_new.u(j);
                walker.theta(j)=walker_new.theta(j);
                reject2 = reject2 + 1 ;
            end
 
       end
         
      
      Gamma= walker.u(lent-1)*exp(-1/walker.theta(1)) ;
      
      walker.u(lent)= Gamma + (( walker.u(lent)-Gamma)*sqrt(1-(Sm*Sm))) + (Sm*sqrt(sigsq1)*randn);
      walker.theta(lent) = walker.u(lent) * walker.u(lent) ;
       F=obs(lent-1)-obs(lent-2);
      walker.logl=walker.logl-log(sqrt(4*pi*walker.theta(lent))) - ((F*F)/(4 * walker.theta(lent)))+log(sqrt(4*pi*walker_new.theta(lent))) + ((F*F)/(4 * walker_new.theta(lent)));  
            
            if(walker.logl > logLstar)
                
                accept2 = accept2 + 1 ;
            else
		walker.logl= walker.logl+log(sqrt(4*pi*walker.theta(lent))) +((F*F)/(4 * walker.theta(lent)))-log(sqrt(4*pi*walker_new.theta(lent))) -((F*F)/(4 * walker_new.theta(lent))); 
                walker.u(lent)=walker_new.u(lent);
                walker.theta(lent)=walker_new.theta(lent);
                reject2 = reject2 + 1 ;
            end
            
        walker_new = walker ;
        
        i = i + 1;
   end
   R1 = reject1 / model.options.nsteps;  % Ratio of rejectance
   step_mod = min(step_mod * exp(0.5-R1),1);	% Update step modifier for all but Y
   R2 = reject2 / (accept2 + reject2);
   stp_mod = min(stp_mod * exp(0.5-R2),1); %update step modifier for Y
   
end

     
