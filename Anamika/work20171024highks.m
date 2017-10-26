%% Run wrapper for run numeric sims for high ks
%changing lambda_hat through Kg
%changing c0_hat = c0/Km
%limit of high ks 
%with nmito = 70 for comparison set

% set up parameter values different from default

l_llim = -1.75;
l_ulim = 1;
nl = 101;
c0_llim = -3;
c0_ulim = 2;
nc0 = 102;
options.gpts = 100;
options.nmito = 70;
options.L = 500;
options.msize = 1;
options.D = 140;
options.dodisplay = 0;
options.dttol = 1e-2;
options.nsteps = 1e4;



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
end
%save the workspace
formatOut = 'yyyymmdd';
date = datestr(datetime('today'),formatOut);
%save workspace with today's date'
filename = strcat('workspace_',date,'gstatnmito70');
save (filename);