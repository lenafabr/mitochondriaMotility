% look at anamika's discrete sim results
load('workspace_20170901edges.mat')
mean(varmito*6/opt.L^2 - 0.5)

load('workspace_20170902continuous_sol.mat')
mean(varmito*6/opt.L^2 - 0.5)

load('workspace_20170901uniform.mat')
mean(varmito*6/opt.L^2 - 0.5)

optdisc = opt;

%% results with 500 mitochondria, same overall glucose consumption rate
load('../results/discretesims_500mito_unif_ks5_20170928.mat')
varmetric = varmito*6/(opt.L-1)^2 - 0.5
mean(varmetric)
optdisc = opt;
%% iterative sims for same parameters
options.ks = optdisc.ks*10; 
options.kg = optdisc.kg*0.1;
options.c0 = optdisc.c0*0.1;
%options.ks = optdisc.ks; 
%options.kg = optdisc.kg;
%options.c0 = optdisc.c0;
options.cend = options.c0;
options.nmito = 500;%optdisc.nmito;

options.dodisplay=1;
options.showevery=100;
options.delt = 0.001;
options.L = 500;

options.gpts = 500;
options.nstep = 1e5;
options.dttol = 1e-3;
[gluc,Tmito,Smito,Smito_int,normdtg,gluc_init,opt,xpos,lmdh,ftc] = runiterativesims(options);

varmito = var(xpos,Tmito);
varmetric = varmito*6/options.L^2 - 0.5

%% run more discrete sims for comparison
options.nstep = 2e5;
options.kw = 1;
options.showmito = 1;
options.showevery = 500;

% TRY STARTING UNIFORM
options.startpos = -1;

options.delt = 0.05;
nitr = 100;
options.pstartwalk = 1;
options.startgluc = gluc;
%
clear varmito gluc_dis mitopos_dis
%
for j = 1:1:nitr
   

    [gluc, mitopos, mitostate, opt] = runmitosim_michaelis2(options);
    varmito(j) = var(mitopos) ; %variance in mitochondria position distribution;
    gluc_dis(j,:) = gluc;
    mitopos_dis(j,:) = mitopos;
    
    [j varmito(j)]
    
    save('/home/ekoslover/proj/mitochondriaMotility/results/discretesims_500mito_unif_cpt1_20170928.mat')
end
