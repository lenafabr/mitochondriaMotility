%% Run wrapper for run numeric sims for high ks
%fixed lambda_hat, c0
%5x space discretization
%limit of high ks 

% set up parameter values different from default

options.L = 500;
options.msize = 1;
options.D = 140;
options.dodisplay = 1;
options.showevery=50000;
options.nmito = 100;
options.dttol = 1e-2;
options.nstep = 1e5;
options.delt=1e-3;
lambda_hat = 0.1;
options.gpts = 100; %5 times 100
options.c0 = 5;
options.kg = options.D ./ (options.nmito * options.msize * options.L * (lambda_hat)^2);
options.cend = options.c0;
[Gstat,Tmito,ksx_stat,gluc_init,opt,xpos,ftc] = gstatsim(options);
var_mito= var(xpos,Tmito) ; %variance in mitochondria position distribution;
varmetric = 6*var_mito/options.L^2 - 0.5
%save the workspace
formatOut = 'yyyymmdd';
date = datestr(datetime('today'),formatOut);
%save workspace with today's date'
%filename = strcat('workspace_',date,'gstat5x');
%save (filename);
