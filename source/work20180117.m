options = struct();
options.ks=500;
options.kg = 1;
options.Km = 1;
options.c0 = 1;
options.cend = 1;
options.ks=5;
options.nmito=100;
options.delt=0.001;
options.showevery=1000;
options.gpts = 100;
[gluc, mitopos, mitostate, opt] = rundiscretesims(options)

%% attempt to reproduce old results using old parameters

options.nstep = 1e5;
options.Km = 0.1;
options.c0 = 0.1;
options.kw = 1;
options.L = 500;
options.ks=100;
options.kg=1;
options.dodisplay=1;
options.showevery=100;
options.nmito=75;

options.delt = 0.05;
options.gpts = 100;
nitr = 100;
%
clear varmito gluc_dis mitopos_dis
for j = 1:1:nitr
    [gluc, mitopos, mitostate, opt] = rundiscretesims(options);
    varmito(j) = var(mitopos) ; %variance in mitochondria position distribution;
    gluc_dis(j,:) = gluc;
    mitopos_dis(j,:) = mitopos;
    
    [j varmito(j)]
end


%% attempt to reproduce old results using old parameters


options.nstep = 1e6;
options.Km = 1;
options.c0 = 1;
options.kw = 1;
options.showmito = 1;
options.showevery = 100;
options.msize = 1;
options.L = 500;
options.ks=5;
options.kg=0.2;
options.dodisplay=1;
options.showevery=500;
options.nmito=100;

% TRY STARTING UNIFORM
options.startpos = -1;
% TRY STARTING AT EDGES
%options.startpos = [10,490];
% TRY STARTING AT CONTINUOUS SOLUTION
%options.startpos = startpos';

options.delt = 0.05;
options.gpts = 100;
nitr = 100;
options.pstartwalk = 1;
%options.startgluc = gluc;
%
clear varmito gluc_dis mitopos_dis
for j = 1:1:nitr
    [gluc, mitopos, mitostate, opt] = rundiscretesims(options);
    varmito(j) = var(mitopos) ; %variance in mitochondria position distribution;
    gluc_dis(j,:) = gluc;
    mitopos_dis(j,:) = mitopos;
    
    [j varmito(j)]
end