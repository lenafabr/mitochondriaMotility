% code to simulate constant permeability case
%constant permeability appears as a sink term in the diffusion equation
function [gluc,Tmito,Smito,Smito_int,normdtg,gluc_init,opt,xpos,lmdh,ftc] = constantpsims(options)
%% set up default simulation parameters
opt = struct();
opt.msize = 1; % mitochondria size

opt.L = 500; % domain size
opt.vel = 1; % mitochondria velocity
opt.D = 140;% glucose diffusion coefficient
opt.kw = 1; % rate of starting a walk
opt.ks = 1; % rate of stopping is ks*[gluc]
opt.kg = 0.1; % rate of glucose consumption (if linear) 
opt.Km = 1;

% starting glucose distribution
% default is to start linear
opt.startgluc = [];

opt.nmito = 75; % number of mitochondria
opt.gpts = 100; % number of discrete spatial points for evaluating gluc concentration
opt.delt = 1e-2; % time-step
opt.nstep = 1e5; % number of steps to run

%External glucose profile parameters
%default is linear external glucose profile
opt.c0 = 1;
opt.cend = 0.1; 
%Permeability term
opt.P = 0.1;
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
% set up dimensionless parameters
%Nondimensionalize by L, L^2/D and Km
tscale = (opt.L)^2 / opt.D;
lscale = opt.L;
cscale = opt.Km;
Lh = opt.L/lscale;
kwh = opt.kw*tscale;
ksh = opt.ks*cscale*tscale; %verify
Dh = opt.D*tscale/lscale^2;
kgh = opt.kg*tscale;
Kmh = opt.Km/cscale;
c0h = opt.c0/cscale;
cendh = opt.cend/cscale;
Ph = opt.P * tscale;
% spatial resolution
dx = Lh/(opt.gpts - 1);

%% Initialize start glucose concentration with analytical solution
%Analytical solution obtained by assuming uniform distribution of
%mitochondria - should I use some other assumption?

% spatial positions at which glucose is evaluated
% index 1 = point on domain edge
xpos = linspace(0,Lh,opt.gpts)';
lmdh = sqrt(Dh./(kgh*opt.nmito*Lh)); %lambda-hat
gluc_init = c0h *  cosh((xpos-Lh/2)./(Lh * lmdh)) ./ cosh(0.5/lmdh);
gluc = gluc_init;
d2g = zeros(opt.gpts,1);
dtg = zeros(opt.gpts,1);
%define external glucose profile (nondimensionalized)
C_out = ((cendh - c0h) * (xpos - Lh/2) / Lh) + (cendh + c0h)/2;

ftc = 0; %flag for failing to converge. Is 1 when fails to converge. 
normdtg = inf;
dtcutoff = opt.dttol;
initglucint = dx * trapz(gluc_init);

%% Iterative process
%continues till steady state
%steady state condition set by time derivative being small enough
step = 0;
while (normdtg > dtcutoff)
    %Calculate distribution of total number of mitochondria
    ksx = ksh * Kmh * gluc ./ (Kmh + gluc);
    ksx_int = dx * trapz(ksx);
    Tmito = (ksx/kwh + 1) ./ (Lh + (ksx_int/kwh));
    Smito = (ksx/kwh) ./ (Lh + (ksx_int/kwh));
    Smito_int = dx * trapz(Smito);
    
    %Calculate the change in glucose concentration
    d2g(2:end-1) = (gluc(3:end)+gluc(1:end-2) - 2*gluc(2:end-1))/dx^2; %space double derivative
    % time derivative of glucose
    dtg(2:end-1) = Dh*d2g(2:end-1) - (kgh * Kmh * opt.nmito * opt.msize) * (gluc(2:end-1) .* Tmito(2:end-1)) ./ (Kmh + gluc(2:end-1)) + (Ph*(C_out(2:end-1) - gluc(2:end-1)));
    normdtg = norm(dtg);
    gluc = gluc+dtg*opt.delt;
    %implement reflecting boundary condition - is this implementation
    %correct?
    gluc(1) = gluc(2);
    gluc(end) = gluc(end-1);
       
    if (any(gluc < -1e-3))
         disp('Concentration went negative. Try smaller timestep.')
        ftc = 1;
        return
    end
    step = step+1;
    
    if (step>opt.nstep)
        disp('Failed to converge')
        ftc = 1;
        return
    end
    
    if (opt.dodisplay && mod(step,opt.showevery)==0)
        
        plot(xpos,gluc_init,'k--')
        hold all
        plot(xpos,gluc,'b.-')
        plot(xpos,C_out,'g--')
        plot(xpos,Tmito*initglucint,'r.-')
        hold off
        drawnow
    end
end
end
 