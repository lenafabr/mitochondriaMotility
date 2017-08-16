options = struct();
options.nmito = 14*5;

%c0 = 0.01

options.L = 500;
options.D = 140;
options.Km = 1;
%options.startpos = 50;
options.pstartwalk = 1;
options.nstep = 1e5;
options.restart = 1;
options.msize = 1;

options.kw = 1;
options.delt=5e-3;

options.dodisplay=1;
options.showevery=1;

A2_llim = -2;
A2_ulim = 2;
nA2 = 5;
l_llim = -2;
l_ulim = 0.5;
nl = 101;
c0_llim = -1;
c0_ulim = 2.5;
nc0 = 102;
%run the function
%change A2 at every iteration, thereby changing ks
%change lambda hat at every iteration, thereby changing kg
%change c0 at each iteration, thereby changing A

c0list = logspace(c0_llim,c0_ulim,nc0);

gluc_all = zeros(100,nl,nc0,nA2);
Smito_all = zeros(100,nl,nc0,nA2);
Smito_int_all = zeros(nl,nc0,nA2);

A2 = 50;
lambda_hat = 1e-2;
options.c0 = 1;
options.ks = A2 * options.kw ./ options.Km;

options.kg = options.D ./ (options.nmito * options.msize * options.L * (lambda_hat^2));

[gluc,Tmito,Smito,Smito_int,normdtg,gluc_init,opt,xpos,lmdh,ftc] = runiterativesims(options);

var_mito = var(xpos,Tmito) ; %variance in mitochondria position distribution;
