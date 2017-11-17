% code to simulate constant permeability case
%constant permeability appears as a sink term in the diffusion equation
function [gluc,Tmito,Smito,Smito_int,normdtg,gluc_init,opt,xpos,lmdh,ftc] = constantpsims(options)
%% set up default simulation parameters
opt = struct();
opt.msize = 1; % mitochondria size

opt.L = 500; % domain size
opt.D = 140;% glucose diffusion coefficient
opt.kw = 1; % rate of starting a walk
opt.ks = 1; % rate of stopping is ks*[gluc]
opt.kg = 0.1; % rate of glucose consumption (if linear) 
opt.Km = 1;

% starting glucose distribution
% default is to start linear
opt.startgluc = [];

opt.nmito = 70; % number of mitochondria
opt.gpts = 100; % number of discrete spatial points for evaluating gluc concentration
opt.delt = 1e-2; % time-step
opt.nstep = 1e5; % number of steps to run

%External glucose profile parameters
%default is linear external glucose profile
opt.c0 = 1;
opt.cend = 0.1; 

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

%Permeability term
opt.P = 0.1;

% copy over supplied options to replace the defaults
if (exist('options')==1)
    opt = copyStruct(options, opt);
end

opt.f = opt.nmito * opt.msize / opt.L;

% set up dimensionless parameters

Lh = opt.L/opt.msize;
kwh = 1;
ksh = opt.ks/opt.kw*opt.c0;
Dh = opt.D/opt.msize^2/opt.kw;
kgh = opt.kg/opt.kw;
Kmh = opt.Km/opt.c0;

% spatial resolution
dx = Lh/(opt.gpts - 1);

%% Initialize start glucose concentration with analytical solution
%Analytical solution obtained by assuming uniform distribution of
%mitochondria 

xpos = linspace(0,Lh,opt.gpts)';


%define external glucose profile, nondimensionalized by c0
C_out = ((opt.cend - opt.c0) * (xpos - Lh/2) / Lh) + (opt.cend + opt.c0)/2;
Cx = C_out / opt.c0;

% spatial positions at which glucose is evaluated
% index 1 = point on domain edge
lmdh = sqrt(Dh./(kgh*opt.nmito*Lh+opt.P*Lh^2)); %lambda-hat
%gluc_init =  cosh((xpos-Lh/2)./(Lh * lmdh)) ./ cosh(0.5/lmdh);
keff = Dh/(lmdh*Lh)^2;
gluc_init = opt.P/keff*Cx - ...
    opt.P*(opt.cend - opt.c0)/opt.c0 * lmdh/keff/cosh(1/lmdh/2)*sinh((xpos-Lh/2)/Lh/lmdh);
gluc = gluc_init;
d2g = zeros(opt.gpts,1);
dtg = zeros(opt.gpts,1);


ftc = 0; %flag for failing to converge. Is 1 when fails to converge. 
normdtg = inf;
dtcutoff = opt.dttol*(kgh*Kmh/(Kmh+1));
spacing = Lh/opt.gpts; %integration spacing
initglucint = spacing * trapz(gluc_init);

%% Iterative process
%continues till steady state
%steady state condition set by time derivative being small enough
step = 0;
while (normdtg > dtcutoff)
    
    %Calculate distribution of total number of mitochondria
    ksx = ksh * Kmh * gluc ./ (Kmh + gluc);
    ksx_int = spacing * trapz(ksx);
    Tmito = (ksx/kwh + 1) ./ (Lh + (ksx_int/kwh));
    Smito = (ksx/kwh) ./ (Lh + (ksx_int/kwh));
    Smito_int = spacing * trapz(Smito);
    
    %Calculate the change in glucose concentration
    d2g(2:end-1) = (gluc(3:end)+gluc(1:end-2) - 2*gluc(2:end-1))/dx^2; %space double derivative
    % time derivative of glucose
    dtg(2:end-1) = Dh*d2g(2:end-1) - (kgh * Kmh * opt.nmito * opt.msize) * (gluc(2:end-1) .* Tmito(2:end-1)) ./ (Kmh + gluc(2:end-1)) + (opt.P*(Cx(2:end-1) - gluc(2:end-1)));
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
        plot(xpos,Cx,'g--')
        %plot(xpos,Tmito*initglucint,'r.-')
        plot(xpos,Tmito*100,'r.-')
        title(sprintf('Step %d; ndtg %f',step,normdtg))        
        hold off
        drawnow
    end
end

% reset glucose concentrations back to real units
gluc = gluc*opt.c0;
gluc_init = gluc_init*opt.c0;
end
 