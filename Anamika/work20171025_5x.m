%% Run wrapper for run numeric sims for high ks
%fixed lambda_hat, c0
%5x space discretization
%limit of high ks 

% set up parameter values different from default

options.L = 500;
options.msize = 1;
options.D = 140;
options.dodisplay = 0;
options.showevery=50000;
options.nmito = 100;
options.dttol = 1e-2;
options.nstep = 1e6;
options.delt=1e-3;
lambda_hat = 0.1;
options.gpts = 500; %5 times 100
%options.c0 = 5;
c0list = 5:5:100;
options.kg = options.D ./ (options.nmito * options.msize * options.L * (lambda_hat)^2);
for i = 1:1:6
    options.c0 = c0list(i);
    options.cend = options.c0;
    [Gstat,Tmito,ksx_stat,gluc_init,opt,xpos,ftc] = gstatsim(options);
    ftc_all(i) = ftc;
    Gstat_all(i,:) = Gstat;
    Tmito_all(i,:) = Tmito;
    var_mito_all(i)= var(xpos,Tmito) ; %variance in mitochondria position distribution;
    varmetric(i) = 6*var_mito_all(i)/options.L^2 - 0.5;
end

%save the workspace
formatOut = 'yyyymmdd';
date = datestr(datetime('today'),formatOut);
%save workspace with today's date'
filename = strcat('workspace_',date,'gstat5x1e6');
save (filename);
