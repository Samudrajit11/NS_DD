function FB_DD_BM_skel()

addpath('util')
data_id2 = '/home/sthapa/NS_FBM/Check/Ensemble_smart_length/delT1/Datafiles/';
misc.data_id= '/home/sthapa/NS_FBM/Check/Ensemble_smart_length/delT1/Resultfiles/';
file_num=10;
for jj = 1:file_num
 
    filename1 = sprintf('%s_%d%s', 'DD_tau5_traj',jj,'.txt');
    data_path = [data_id2,filename1];
%data_path = [misc.data_id,'subdiffusive_track2.txt'];

%loading data 
    data = load(data_path);
    tau_length=length(data);
    data = transpose(data);

%Subtract points to obtain steps
    obs = data(:,2:end) - data(:,1:end-1);

%logl_true = fbm_logl(obs,[20 0 0 10 0.75],[1 1 1 1])



%Specify framerate and pixelsize
%p_size = 1;  %Pixelsize in data in micrometer/pixel
    tau = 1;     %Time between data point in seconds

%Specify prior ranges
    sigmaHmin=10^(-3);  %Minimum value for step-deviation in micrometer^2
    sigmaHmax=10^(3);      %Maximum value for step deviation in micrometer^2

    vmax=10;      %Maximum absolute value of drift in micrometer/s

    noisemin=0;  %Minimum value for noise parameter in micrometers
    noisemax=1;  %Maximum value for noise parameter in micrometers

    Hmin=0;      %Minimum value for Hurst parameter
    Hmax=1;      %Maximum value for Hurst parameter

%Convert prior ranges to dimensions of pixels and frames for convenience 
% ranges - a 2x5 array of minimum and maximum values for the 5 parameter:
%   - diffusion coefficient, drift velocity (in two dimentions), noise
%     parameter, and Hurst exponent.
    ranges2=[sigmaHmin sigmaHmax;...
        -vmax vmax ;-vmax vmax;...
        noisemin noisemax; Hmin Hmax];

    N=length(data(1,:));
    len=N+1;
%Specify prior ranges
%Dmin=10^(-3);
%Dmax=10^(3);   

    tmin=1;
    tmax=tau_length*tau;   

    dsmin = 10^(-2);
    dsmax = 10^(1);

%Convert prior ranges to dimensions of pixels and frames for convenience 
% ranges - a 2x5 array of minimum and maximum values for the 4 parameters:
%   - Mobility coefficient, mobility exponent, force, and measurement error
    ranges=[tmin tmax ; dsmin dsmax ];


%Specify options
    options.nwalkers=200; % Number of walkers to be generated
    options.stoprat=10^(-4);
    options.nsteps=50;
%    options.maxsamples=1000;
%options.nlist= [1 2 4 8 16 32];
%options.trackmax = 100;


%Specify the models
    final=@(x) x(end); % extract last element in a vector
    MM = [0 0 0 0; 0 0 1 0; 0 0 0 1; 0 0 1 1 ];
    for i=1:4
        models(i).genu=@() util_generate_u(sum(MM(i,:))+1);
        models(i).options=options;
        models(i).logl=@(obs,theta) fbm_logl(obs,fbm_params(theta,MM(i,:)),MM(i,:));
        models(i).invprior=@(u) fbm_invprior(u,ranges2,MM(i,:));
        %models(i).scaling = @(obs,n) fbm_scaling(obs,n);
        models(i).evolver=@(obs,model,logLstar,walker,step_mod)  ns_evolve_exp(obs,model,logLstar,walker,step_mod);
        %models(i).replicate = @(obs,theta,n) fbm_replicate(obs,fbm_params(theta,MM(i,:)),n);
        %models(i).logl_n = @(obs,theta,n) fbm_logl_n(obs,fbm_params(theta,MM(i,:)),MM(i,:),n);
        models(i).labels=[1];
        for j=1:length(MM(i,:))
            if MM(i,j)==1
            models(i).labels=[models(i).labels j+1];
            end
        end
        models(i).add{1}=@(theta) theta(1)^2/(2*tau^(2*final(fbm_params(theta,MM(i,:)))));
        models(i).labels=[models(i).labels 6];
        for j=1:2
            if MM(i,j)==1
            models(i).add{end+1}=@(theta) theta(j+1)/tau;
            models(i).labels=[models(i).labels 6+j];
            end
        end
    end

    models(5).options=options;
    models(5).genu=@() geny(ranges, len);
    models(5).logl=@(obs,theta) DD_logl_m(obs,theta);
    models(5).invprior=@(u) DDU_invprior(u,ranges);
    %models.u_evolve=@(u, delta_u)  u_evolve(u,delta_u,ranges);
    %models.evolver=@ns_evolve_y;
    models(5).evolver=@(obs,model,logLstar,walker,step_mod,stp_mod)  ns_new_evolve_y_exp(obs,model,logLstar,walker,step_mod,stp_mod,ranges);
    models(5).labels=[9:10 6*ones(1,(N-1))];

%Labels for the parameters
    misc.labels=...
    ['Step deviation: ';...
    'x-bias/step:    ';...
    'y-bias/step:    ';...
    'Measurem. err.: ';...
    'Hurst exponent: ';...
    'D_H constant:   ';...
    'x-velocity:     ';...
    'y-velocity:     ';...
    'Tau:            ';...
    'D_star:         '];

%Tell ns_processdataset to write a summary file
    filename2 = sprintf('%s_%d%s','nomaxResults_delT1_DD_tau5_traj',jj,'.txt' );
    misc.nssummary=filename2;
    
    filename3=sprintf('%s_%d','nomaxResults_delT1_DD_tau5_traj',jj);
    path=fullfile(misc.data_id,filename3);

    [results] = ns_processdataset(obs,models,misc);
    save(path,'results')
    save(path,'data','-append')
    save(path,'options','-append')
    save(path,'ranges','-append')
end
end

