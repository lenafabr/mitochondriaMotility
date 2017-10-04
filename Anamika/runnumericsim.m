function [Gstat,Gstat_prime,Tmito,ksx_calc,gluc_init,opt,xpos] = runnumericsim(options)
%function to get a stationary soltuion for G(x) in the limit of high Ks
%Numerically solve the differential equation to get G(x)
%Then get T(x) from this G_stat(x)
%% set up default simulation parameters
opt = struct();

opt.kg = 1; % rate of glucose consumption
opt.c0 = 0.1; % fixed glucose concentration
opt.msize = 1; % mitochondria size

opt.L = 500; % domain size
opt.vel = 1; % mitochondria velocity
opt.D = 140;% glucose diffusion coefficient
opt.kw = 1; % rate of starting a walk
opt.ks = 50; % rate of stopping is ks*[gluc]
opt.kg = 0.1; % rate of glucose consumption (if linear)
opt.Km = 1;

% starting glucose distribution
% default is to start linear
opt.startgluc = [];

opt.nmito = 100; % number of mitochondria
opt.gpts = 100; % number of discrete spatial points for evaluating gluc concentration
opt.delt = 1e-4; % time-step
opt.nstep = 1e5; % number of steps to run



% boundary conditions on the far size
% positive = fixed concentration at the boundary
% negative = reflecting boundary
opt.cend = opt.c0;

% starting position of mitochondria
% default (<0) means start uniformly
opt.startpos = -1;
% prob mitochondria start in walking state
% default is use equilibrium probability at c0 conc
opt.pstartwalk = opt.kw/(opt.kw + opt.ks*opt.c0);

% tolerance for "small time derivative"
opt.dttol = 1e-4;

% displaying plots
opt.dodisplay = 1;
opt.showevery = 1;

opt.restart = 1; % flag to enable continuing previous sims

% copy over supplied options to replace the defaults
if (exist('options')==1)
    opt = copyStruct(options, opt);
end

opt.f = opt.nmito * opt.msize / opt.L;

% set up dimensionless parameters

Lh = opt.L/opt.msize;
% velh = opt.vel/opt.msize/opt.kg;
% kwh = opt.kw/opt.kg;
% ksh = opt.ks/opt.kg*opt.c0;
% Dh = opt.D/opt.msize^2/opt.kg;
velh = 1;
kwh = opt.kw/opt.vel*opt.msize;
ksh = opt.ks/opt.vel*opt.c0*opt.msize;
Dh = opt.D/opt.msize/opt.vel;
kgh = opt.kg*opt.msize/opt.vel;
Kmh = opt.Km/opt.c0;

% spatial resolution
dx = Lh/(opt.gpts - 1);

%% Initialize start glucose concentration with analytical solution
%Analytical solution obtained by assuming uniform distribution of
%mitochondria

% spatial positions at which glucose is evaluated
% index 1 = point on domain edge
xpos = linspace(0,Lh,opt.gpts)';
lmdh = sqrt(Dh./(kgh*opt.nmito*Lh)); %lambda-hat
gluc_init =  cosh((xpos-Lh/2)./(Lh * lmdh)) ./ cosh(0.5/lmdh);
gluc = gluc_init;
d2g = zeros(opt.gpts,1);
% dtg = zeros(opt.gpts,1);
% % d2g(2:end-1) = (gluc(3:end)+gluc(1:end-2) - 2*gluc(2:end-1))/dx^2; %space double derivative
% %dtg(1) = 0;
% %dtg(end) = 0; %fixed boundary conditions
% % dtg(2:end-1) = Dh*d2g(2:end-1);
% % dtg_init = dtg;
% ftc = 0; %flag for failing to converge. Is 1 when fails to converge.
%
% normdtg = inf;

dtcutoff = opt.dttol*(kgh*Kmh/(Kmh+1));
spacing = Lh/opt.gpts; %integration spacing

initglucint = spacing * trapz(gluc_init);


%% Solve the DE for G using bvp - using bvp4c
%G_prime_init = (gluc_init(2) - gluc_init(1)) / spacing;
a = gluc_init(1);
b = gluc_init(2);
ksx = ksh * Kmh * gluc ./ (Kmh + gluc);
ksx_int = spacing * trapz(ksx);
const = opt.nmito * opt.msize ./ (Dh * ksx_int);
dGdx = @(x,G) [G(2); const * (((kgh*Kmh*G(1))/ (Kmh + G(1)))^2)];
res = @(G0,GL) [G0(1)-a; GL(1)-b];
solinit = bvpinit(xpos,[opt.c0,0]);
Gstat_fn = bvp4c(dGdx,res,solinit);
Gstat_pts = deval(xpos,Gstat_fn);
Gstat = Gstat_pts(1,:);
Gstat_prime = Gstat_pts(2,:);
ksx_calc = ksh * Kmh * Gstat ./ (Kmh + Gstat);
ksx_calc_int = spacing * trapz(ksx_calc);
Tmito = ksx_calc ./ ksx_calc_int;


%% Solve the DE w/o using inbuilt function
% so that we can set our own tolerance values

end


