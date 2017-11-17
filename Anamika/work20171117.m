%comparing constant permeability sims in the linear limit
%% Define values of parameters
options = struct();
options.msize = 1; % mitochondria size
options.L = 500; % domain size
options.vel = 1; % mitochondria velocity
options.D = 140;% glucose diffusion coefficient
options.kw = 1; % rate of starting a walk
options.ks = 1; % rate of stopping is ks*[gluc]
options.kg = 100; % rate of glucose consumption (if linear) 
options.Km = 10; %linear regime
options.nmito = 70; % number of mitochondria
options.gpts = 100; % number of discrete spatial points for evaluating gluc concentration
options.delt = 1e-2; % time-step
options.nstep = 1e5; % number of steps to run

options.delt = 1e-2; % time-step
options.nstep = 1e5; % number of steps to run
%External glucose profile parameters
%default is linear external glucose profile
options.c0 = 1;
options.cend = 0.1; 

% starting position of mitochondria
% default (<0) means start uniformly
options.startpos = -1;
% prob mitochondria start in walking state
% default is use equilibrium probability at c0 conc
options.pstartwalk = options.kw/(options.kw + options.ks*options.c0);

% tolerance for "small time derivative"
options.dttol = 1e-4;

% displaying plots
options.dodisplay = 0;
options.showevery = 100;

options.restart = 1; % flag to enable continuing previous sims

%Permeability term
options.P = 0.1;

% copy over supplied options to replace the defaults
if (exist('options')==1)
    options = copyStruct(options, options);
end

options.f = options.nmito * options.msize / options.L;

% set up dimensionless parameters

Lh = options.L/options.msize;
velh = 1;
kwh = options.kw/options.vel*options.msize;
ksh = options.ks/options.vel*options.c0*options.msize;
Dh = options.D/options.msize/options.vel;
kgh = options.kg*options.msize/options.vel;
Kmh = options.Km/options.c0;
Ph = options.P * options.msize / options.vel;

xpos = linspace(0,Lh,options.gpts)';

% spatial resolution
dx = Lh/(options.gpts - 1);

%% First compare simulation with analytical results
[gluc,Tmito,Smito,Smito_int,normdtg,gluc_init,opt,xpos,lmdh,ftc] = constantpsims(options);

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

