%run wrapper for runmitosim_michaelis2 function - with changing nmito
%Discrete case - compare with the continuous case

%define options structure

options = struct();

options.L = 500;
options.D = 140;
%options.Km = 0.1/c0;
%options.startpos = 50;
options.pstartwalk = 1;
options.nstep = 1e6;
options.restart = 1;
options.km = 1;
options.c0 = 1;

%options.ks = (1/4.8*1e-6)*c0*(10^-3*6e23/1000/1e12*4^2);
options.kw = 1 * 0.01;
options.delt=5e-2;
options.msize = 1; % mitochondria size
options.dodisplay=0;
options.showevery=1;
options.ks = 1;
options.kg = 0.2;


%run the function
%change nmito in every iteration and calculate the variance

for i = 1:1:1000
    options.nmito = i;
    [gluc, mitopos, mitostate, opt,lmdh] = runmitosim_michaelis2(options)
    var_metric(i) = var(mitopos) ; %variance in mitochondria position distribution;
end