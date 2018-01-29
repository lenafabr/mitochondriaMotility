%% figure out iterative vs discrete comparison: are we converging on the discrete sims?

options.nstep = 1e6;
options.Km = 0.1;
options.c0 = 0.1;
options.kw = 1;
options.L = 500;
options.ks=100;
options.kg=1/10;
options.dodisplay=0;
options.showevery=100;
options.nmito=750;

options.delt = 1e-5;
options.nstep = 1e7;
options.gpts = 100;
[gluc_itr,Tmito_itr,Smito_itr,Smito_int_itr,normdtg,gluc_init,opt,xpos,lmdh,ftc] = runiterativesims(options);



%% Run discrete simulations starting with iterative solution (for comparison with analytical results)
options=struct();
options.nstep = 2e5;
options.Km = 0.1;
options.c0 = 0.1;
options.kw = 1;
options.L = 500;
options.ks=100;
options.kg=1;
options.dodisplay=0;
options.showevery=500;
options.nmito=75;
options.startgluc = gluc_itr;

options.delt = 0.05;
options.gpts = 100;
nitr = 100;

clear varmito gluc_dis mitopos_dis
for j = 1:1:nitr
    [gluc, mitopos, mitostate, opt] = rundiscretesims(options);
    varmito(j) = var(mitopos) ; %variance in mitochondria position distribution;
    gluc_dis(j,:) = gluc;
    mitopos_dis(j,:) = mitopos;
    
    [j varmito(j)]
    
    %save the workspace
    formatOut = 'yyyymmdd';
    date = datestr(datetime('today'),formatOut);
    %save workspace with today's date'
    filename = strcat('/home/ekoslover/proj/mitochondriaMotility/results/workspace_',date,'discretesims_100itr_startgluciter');
    save (filename);
end