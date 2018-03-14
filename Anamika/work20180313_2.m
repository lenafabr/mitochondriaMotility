%For uniform permeability, using a conc profile b/w 50Km and 0.1 Km
%Surface plot of avg mito position and enrichment around high conc point 
%vs lambda_hat and P

%% Define values of parameters
options = struct();
options.ks = 1e3; % rate of stopping is ks*[gluc], high ks regime
options.Km = 0.1; 
options.D = 140;
options.nmito = 75; % number of mitochondria
options.msize = 1;
options.L = 500;
options.gpts = 100; % number of discrete spatial points for evaluating gluc concentration
options.delt = 1e-3; % time-step
options.nstep = 1e6; % number of steps to run
%default is linear external glucose profile
options.c0 = 50; %Using exptl paper
options.cend = 0.1; %Exptl
% tolerance for "small time derivative"
options.dttol = 1e-2;
options.nstep = 1e8;
options.delt = 1e-6;
% displaying plots
options.dodisplay = 0;
options.showevery = 100;

%% run wrapper with varying P and lambda_hat
lh_ulim = -2;
lh_llim = 2;
P_ulim = 3;
P_llim = -3;
nlh = 11;
nP = 12;
lhlist = logspace(lh_llim,lh_ulim,nlh);
Plist = logspace(P_llim,P_ulim,nP);
for i = 1:1:nlh
    lambda_hat = lhlist(i);
    options.kg = options.D ./ (options.nmito * options.msize * options.L * (lambda_hat^2));
    for j = 1:1:nP
        options.P = Plist(j);
        [gluc,Tmito,Smito,Smito_int,normdtg,gluc_init,opt,xpos,lmdh,ftc] = constantpsims(options);
        ftc_matrix(i,j) = ftc;
        option_list(i,j) = options;
        gluc_init_all(:,i,j) = gluc_init;
        gluc_all(:,i,j) = gluc;
        Tmito_all(:,i,j) = Tmito;
        Smito_all(:,i,j) = Smito;
        Percent_complete = (i/nlh) * 100 
    end
end
%save the workspace
formatOut = 'yyyymmdd';
date = datestr(datetime('today'),formatOut);
%save workspace with today's date'
filename = strcat('workspace_',date,'p_lh_12_11');
save (filename);

