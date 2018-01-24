%% Run wrapper for run numeric sims for high ks

%changing c0_hat = c0/Km
%changing lambda_hat
%nmito = 75

% set up parameter values different from default
lambda_hat = 0.06;
options.D = 140;
c0_llim = -2;
c0_ulim = 2;
nc0 = 102;
l_llim = -2;
l_ulim = 2;
nl = 101;
options.D = 140;
options.msize = 1;
options.nmito = 75;
options.gpts = 100;
options.L = 500;

options.dodisplay = 0;
options.dttol = 1e-3;
options.delt = 1e-3;
options.nstep = 1e7;


%run the function
%change lambda hat at every iteration, thereby changing kg
%change c0 at each iteration, thereby changing A

c0list = logspace(c0_llim,c0_ulim,nc0);
llist = logspace(l_llim,l_ulim,nl);

Gstat_all = zeros(options.gpts,nl,nc0);
Tmito_all = zeros(options.gpts,nl,nc0);

for i = 1:1:nl
    lambda_hat(i) = llist(i);
    options.kg = options.D ./ (options.nmito * options.msize * options.L * (lambda_hat(i)^2));
    for j = 1:1:nc0
        options.c0 = c0list(j);
        options.cend = options.c0;
        [Gstat,Tmito,ksx_stat,gluc_init,opt,xpos,ftc] = gstatsim(options);
        ftc_matrix(i,j) = ftc;
        option_list(i,j) = opt;
        gluc_init_all(:,i,j) = gluc_init;
        Gstat_all(:,i,j) = Gstat;
        Tmito_all(:,i,j) = Tmito;
        var_mito(i,j) = var(xpos,Tmito) ; %variance in mitochondria position distribution;
        varmetric(i,j) = 6*var_mito(i,j)/options.L^2 - 0.5;
        
    end
    percent_completed_highks = (i/nl) * 100
end
%save the workspace
formatOut = 'yyyymmdd';
date = datestr(datetime('today'),formatOut);
%save workspace with today's date'
filename = strcat('workspace_',date,'gstatnmito75');
save (filename);