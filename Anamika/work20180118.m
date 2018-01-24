%run iterative siumlations for several c0
%for high ks
%This is for the plot showing varying Tmito for c0, blue to red low c0 to
%high c0, with analytical solution shown 
lambda_hat = 0.06;
options.ks = 200;
options.D = 140;
c0_llim = -2;
c0_ulim = 2;
nc0 = 20;
options.D = 140;
options.msize = 1;
options.nmito = 75;
options.gpts = 100;
options.L = 500;
options.kg = 1;
options.Km = 0.1;

options.dodisplay = 0;
options.dttol = 1e-3;
options.delt = 1e-5;
options.nstep = 1e7;

%run the function
%change lambda hat at every iteration, thereby changing kg
%change c0 at each iteration, thereby changing A

c0list = logspace(c0_llim,c0_ulim,nc0);
gluc_all = zeros(options.gpts,nc0);
Tmito_all = zeros(options.gpts,nc0);



for j = 1:1:nc0
    options.c0 = c0list(j);
    options.cend = options.c0;
    [gluc,Tmito,Smito,Smito_int,normdtg,gluc_init,opt,xpos,lmdh,ftc] = runiterativesims(options);
    normdtg_matrix(j) = normdtg;
    ftc_matric(j) = ftc;
    option_list(j) = opt;
    gluc_init_all(:,j) = gluc_init;
    gluc_all(:,j) = gluc;
    Tmito_all(:,j) = Tmito;
    Smito_int_all(:,j) = Smito_int;
    Smito_all(:,j) = Smito;
    var_mito(j) = var(xpos,Tmito) ; %variance in mitochondria position distribution;
    varmetric(j) = 6*var_mito(j) - 0.5;
    percent_completed_c0 = (j/nc0 * 100)
end


%save the workspace
formatOut = 'yyyymmdd';
date = datestr(datetime('today'),formatOut);
%save workspace with today's date'
filename = strcat('workspace_',date,'c0_ks_200');
save (filename);