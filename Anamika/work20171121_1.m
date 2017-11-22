%For uniform permeability, using a conc profile b/w 50Km and 0 Km
%Surface plot of avg mito position and enirchment around high conc point 
%vs Kg and P

%% Define values of parameters
options = struct();
options.ks = 1e3; % rate of stopping is ks*[gluc], high ks regime
options.Km = 1; 
options.nmito = 70; % number of mitochondria
options.gpts = 100; % number of discrete spatial points for evaluating gluc concentration
options.delt = 1e-3; % time-step
options.nstep = 1e6; % number of steps to run
%default is linear external glucose profile
options.c0 = 50; %Using exptl paper
options.cend = 0.1; %Exptl
% tolerance for "small time derivative"
options.dttol = 1e-4;

% displaying plots
options.dodisplay = 0;
options.showevery = 100;

%% run wrapper with varying P and kg
kg_ulim = -1;
kg_llim = 1;
P_ulim = 3;
P_llim = -3;
nkg = 101;
nP = 102;
kglist = logspace(kg_llim,kg_ulim,nkg);
Plist = logspace(P_llim,P_ulim,nP);
for i = 1:1:nkg
    options.kg = kglist(i);
    for j = 1:1:nP
        options.P = Plist(j);
        [gluc,Tmito,Smito,Smito_int,normdtg,gluc_init,opt,xpos,lmdh,ftc] = constantpsims(options);
        ftc_matrix(i,j) = ftc;
        option_list(i,j) = options;
        gluc_init_all(:,i,j) = gluc_init;
        gluc_all(:,i,j) = gluc;
        Tmito_all(:,i,j) = Tmito;
        Smito_all(:,i,j) = Smito;
        Percent_complete = (i/nkg) * 100 
    end
end

