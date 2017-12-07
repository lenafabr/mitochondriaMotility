%% old non-dimensionalization

options.kg = 0.1;
options.c0 = 300;
options.Km = 1;
options.cend = options.c0;
options.dodisplay=1;
options.showevery=1000;
options.delt = 1e-4;
options.ks = 1000;
options.kw=1;
[gluc,Tmito,Smito,Smito_int,normdtg,gluc_init,opt,xpos,lmdh,ftc] = runiterativesims_oldim(options);

var_mito = var(xpos,Tmito) ; %variance in mitochondria position distribution;
varmetric = 6*var_mito/opt.L^2 - 0.5;

%% check new non-dimensionalization for fixed end concentrations
% time units L^2/D, length units L

options.kg = 0.1;
options.c0 = 300;
options.cend = options.c0;
options.dodisplay=1;
options.showevery=1000;
options.delt = 1e-6;
options.ks = 1000;
[gluc,Tmito,Smito,Smito_int,normdtg,gluc_init,opt,xpos,lmdh,ftc] = runiterativesims(options);

var_mito = var(xpos,Tmito) ; %variance in mitochondria position distribution;
varmetric = 6*var_mito - 0.5