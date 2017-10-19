%run wrapper for permeablesims - implementing permeability conditions
%% Run wrapper for changing ks and c0
%lambda_hat fixed at 0.14
%changing c0_hat = c0/Km
%changing ks
%nmito = 100

% set up parameter values different from default
lambda_hat = 0.14;
options.D = 140;
c0_llim = -3;
c0_ulim = 2;
nc0 = 102;
ks_llim = -2;
ks_ulim = 2;
nks = 101;
options.gpts = 100;
options.nmito = 100;
options.L = 500;
options.msize = 1;
options.kg = options.D ./ (options.nmito * options.msize * options.L * (lambda_hat^2));
options.dodisplay = 0;
options.dttol = 1e-3;
options.nstep = 1e4;
options.P = 0.5;

%run the function
%change lambda hat at every iteration, thereby changing kg
%change c0 at each iteration, thereby changing A

c0list = logspace(c0_llim,c0_ulim,nc0);
kslist = logspace(ks_llim,ks_ulim,nks);


for i = 1:1:nks
    options.ks = kslist(i);
    for j = 1:1:nc0
        options.c0 = c0list(j);
        options.cend = options.c0;
        [gluc,Tmito,Smito,Smito_int,normdtg,gluc_init,opt,xpos,lmdh,ftc] = permeablesims(options);
        ftc_matrix(i,j) = ftc;
        option_list(i,j) = opt;
        gluc_init_all(:,i,j) = gluc_init;
        gluc_all(:,i,j) = gluc;
        Tmito_all(:,i,j) = Tmito;
        Smito_int_all(:,i,j) = Smito_int;
        Smito_all(:,i,j) = Smito;
        var_mito(i,j) = var(xpos,Tmito(2:end-1)) ; %variance in mitochondria position distribution;
        varmetric(i,j) = 6*var_mito(i,j)/options.L^2 - 0.5;
    end
    percent_complete = (i/nks) * 100;
    
end


%save the workspace
formatOut = 'yyyymmdd';
date = datestr(datetime('today'),formatOut);
%save workspace with today's date'
filename = strcat('workspace_',date,'perm_ks_c0_05');
save (filename);