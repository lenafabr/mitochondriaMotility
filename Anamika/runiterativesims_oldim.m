function [gluc,Tmito,Smito,Smito_int,normdtg,gluc_init,opt,xpos,lmdh,ftc] = runiterativesims_olddim(options)
%% set up default simulation parameters
opt = struct();

opt.kg = 1; % rate of glucose consumption
opt.c0 = 1; % fixed glucose concentration
opt.msize = 1; % mitochondria size

opt.L = 500; % domain size
opt.vel = 1; % mitochondria velocity
opt.D = 140;% glucose diffusion coefficient
opt.kw = 1; % rate of starting a walk
opt.ks = 1; % rate of stopping is ks*[gluc]
opt.Km = 1;

% starting glucose distribution
% default is to start linear
opt.startgluc = [];

opt.nmito = 1; % number of mitochondria
opt.gpts = 100; % number of discrete spatial points for evaluating gluc concentration
opt.delt = 1e-4; % time-step
opt.nstep = 1e5; % number of steps to run

% boundary conditions on the far size
% positive = fixed concentration at the boundary
% negative = reflecting boundary
opt.cend = 1;
cendh = opt.cend/opt.Km;

% starting position of mitochondria
% default (<0) means start uniformly
opt.startpos = -1;
% tolerance for "small time derivative"
opt.dttol = 1e-4;

% displaying plots
opt.dodisplay = 1;
opt.showevery = 1000;

opt.restart = 1; % flag to enable continuing previous sims

% copy over supplied options to replace the defaults
if (exist('options')==1)
    opt = copyStruct(options, opt);
end
opt.f = opt.nmito * opt.msize / opt.L;
% set up dimensionless parameters
%non-dimensionalize by Km instead of c0
Lh = opt.L/opt.msize
velh = 1;
kwh = opt.kw/opt.vel*opt.msize
ksh = opt.ks/opt.vel*opt.Km*opt.msize
Dh = opt.D/opt.msize/opt.vel
kgh = opt.kg*opt.msize/opt.vel
c0h = opt.c0/opt.Km

% spatial resolution
dx = Lh/(opt.gpts - 1);

%% Initialize start glucose concentration with analytical solution
%Analytical solution obtained by assuming uniform distribution of
%mitochondria

% spatial positions at which glucose is evaluated
% index 1 = point on domain edge
xpos = linspace(0,Lh,opt.gpts)';
lmdh = sqrt(Dh./(kgh*opt.nmito*Lh)); %lambda-hat
gluc_init = c0h * cosh((xpos-Lh/2)./(Lh * lmdh)) ./ cosh(0.5/lmdh);
gluc = gluc_init;
d2g = zeros(opt.gpts,1);
dtg = zeros(opt.gpts,1);
% d2g(2:end-1) = (gluc(3:end)+gluc(1:end-2) - 2*gluc(2:end-1))/dx^2; %space double derivative
%dtg(1) = 0;
%dtg(end) = 0; %fixed boundary conditions
% dtg(2:end-1) = Dh*d2g(2:end-1);
% dtg_init = dtg;
ftc = 0; %flag for failing to converge. Is 1 when fails to converge. 

normdtg = inf;

dtcutoff = opt.dttol*(kgh*opt.Km/(opt.Km+1));
spacing = Lh/opt.gpts; %integration spacing

initglucint = spacing * trapz(gluc_init);

%% Iterative process
%continues till steady state
%steady state condition set by time derivative being small enough
step = 0;
while (normdtg > dtcutoff)
    
    %Calculate distribution of total number of mitochondria
    ksx = ksh * opt.Km * gluc ./ (opt.Km + gluc);
    ksx_int = spacing * trapz(ksx);
    Tmito = (ksx/kwh + 1) ./ (Lh + (ksx_int/kwh));
    Smito = (ksx/kwh) ./ (Lh + (ksx_int/kwh));
    Smito_int = spacing * trapz(Smito);
    
    %Calculate the change in glucose concentration
    d2g(2:end-1) = (gluc(3:end)+gluc(1:end-2) - 2*gluc(2:end-1))/dx^2; %space double derivative
    % time derivative of glucose
    dtg(2:end-1) = Dh*d2g(2:end-1) - (kgh * opt.Km * opt.nmito * opt.msize) * (gluc(2:end-1) .* Tmito(2:end-1)) ./ (opt.Km + gluc(2:end-1));
    normdtg = norm(dtg);
    gluc = gluc+dtg*opt.delt;
       
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
        plot(xpos,Tmito*100,'r.-')
        hold off
        drawnow
    end
end
end

