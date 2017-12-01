%Simulate data for phase diagrams with P = 0.1 
%compare with phase diagrams for fixed end concentrations

%% options for both data sets
l_llim = -1.75;
l_ulim = 1;
nl = 101;
c0_llim = -3;
c0_ulim = 2;
nc0 = 102;
options.gpts = 100;
options.L = 500;
options.msize = 1;
options.D = 140;
options.dodisplay = 0;
options.nmito = 70;
options.dttol = 1e-2;
options.nsteps = 1e6;
options.ks = 20; %kshat = 20
options.kw = 1;
options.Km = 1;
options.P = 0.1;

c0list = logspace(c0_llim,c0_ulim,nc0);
llist = logspace(l_llim,l_ulim,nl);

%% Generate data set for P = 0.1 and only ends permeable
%only need c0 here
for i = 1:1:nl
    lambda_hat(i) = llist(i);
    options.kg = options.D ./ (options.nmito * options.msize * options.L * (lambda_hat(i)^2));
    for j = 1:1:nc0
        options.c0 = c0list(j);
        options.cend = options.c0;
        [gluc,Tmito,Smito,Smito_int,normdtg,gluc_init,opt,xpos,lmdh,ftc] = permeablesims(options);
        ftc_matrix(i,j) = ftc;
        option_list(i,j) = opt;
        gluc_init_all(:,i,j) = gluc_init;
        gluc_all(:,i,j) = gluc;
        Tmito_all(:,i,j) = Tmito;
        Smito_all(:,i,j) = Smito;
        Smito_int_all(i,j) = Smito_int;
        var_mito(i,j) = var(xpos,Tmito(2:end-1)) ; %variance in mitochondria position distribution;
        varmetric(i,j) = 6*var_mito(i,j)/options.L^2 - 0.5;
    end
    percent_complete = (i/nl) * 100  
end

%save the workspace
formatOut = 'yyyymmdd';
date = datestr(datetime('today'),formatOut);
%save workspace with today's date'
filename = strcat('workspace_',date,'p_01_l_c0');
save (filename);

%% Generate data set for fixed concentration at edges

for i = 1:1:nl
    lambda_hat(i) = llist(i);
    options.kg = options.D ./ (options.nmito * options.msize * options.L * (lambda_hat(i)^2));
    for j = 1:1:nc0
        options.c0 = c0list(j);
        options.cend = options.c0;
        [gluc,Tmito,Smito,Smito_int,normdtg,gluc_init,opt,xpos,lmdh,ftc] = runiterativesims(options);
        ftc_matrix(i,j) = ftc;+
        option_list(i,j) = opt;
        gluc_init_all(:,i,j) = gluc_init;
        gluc_all(:,i,j) = gluc;
        Tmito_all(:,i,j) = Tmito;
        Smito_all(:,i,j) = Smito;
        Smito_int_all(i,j) = Smito_int;
        var_mito(i,j) = var(xpos,Tmito) ; %variance in mitochondria position distribution;
        varmetric(i,j) = 6*var_mito(i,j)/options.L^2 - 0.5;
    end
    percent_complete = (i/nl) * 100  
end

%save the workspace
formatOut = 'yyyymmdd';
date = datestr(datetime('today'),formatOut);
%save workspace with today's date'
filename = strcat('workspace_',date,'fixedc0');
save (filename);



