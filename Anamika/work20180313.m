%comparing constant permeability sims in the linear limit
%% Define values of parameters
options = struct();
options.msize = 1; % mitochondria size
options.L = 500; % domain size
options.D = 140;% glucose diffusion coefficient
options.kw = 1; % rate of starting a walk
options.ks = 1e3; % rate of stopping is ks*[gluc]


options.delt = 1e-5; % time-step
options.nstep = 1e7; % number of steps to run
%External glucose profile parameters
%default is linear external glucose profile
options.c0 = 50;
options.cend = 0.1; 


%% First compare simulation with analytical results
[gluc,Tmito,Smito,Smito_int,normdtg,gluc_init,opt,xpos,lmdh,ftc] = constantpsims(options);
Lh = options.L/options.L;
kgh = opt.kg*(options.L^2 /options.D);
Ph = opt.P * (options.L^2 /options.D);
Dh = 1;

kh = kgh * opt.nmito * opt.msize / Lh;
a = Ph / (kh+Ph) * (options.cend-options.c0)/(Lh );
b = Ph / (kh+Ph) * (options.cend+options.c0)/(2);
Beta = sqrt(kh+Ph/Dh);
A = a / (2* cosh(Beta*Lh/2));
gluc_calc = A*exp(-Beta * (xpos-Lh/2)) - A*exp(-Beta * (xpos - Lh/2)) + a*(xpos-Lh/2) + b;

plot(xpos,gluc,'r');
hold on;
plot(xpos,gluc_calc,'b');
hold off

%% Then plug in values from the dataset in the exptl paper
