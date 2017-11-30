function [gluc,Tmito,Smito,Smito_int,normdtg,gluc_init,opt,xpos,lmdh,ftc] = permeablesims(options)
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
opt.delt = 1e-3; % time-step
opt.nstep = 1e5; % number of steps to run

opt.P = 0.1; %permeability
opt.f = opt.nmito * opt.msize / opt.L;
% tolerance for "small time derivative"
opt.dttol = 1e-3;

% displaying plots
opt.dodisplay = 1;
opt.showevery = 1;

opt.restart = 1; % flag to enable continuing previous sims

% copy over supplied options to replace the defaults
if (exist('options')==1)
    opt = copyStruct(options, opt);
end

% set up dimensionless parameters

Lh = opt.L/opt.msize;
velh = 1;
kwh = opt.kw/opt.vel*opt.msize;
ksh = opt.ks/opt.vel*opt.Km*opt.msize;
Dh = opt.D/opt.msize/opt.vel;
kgh = opt.kg*opt.msize/opt.vel;
c0h = opt.c0/opt.Km;
Ph = opt.P * opt.msize / opt.vel;

% spatial resolution
dx = Lh/(opt.gpts - 1);

%% Initialize start glucose concentration with analytical solution
%Analytical solution obtained by assuming uniform distribution of
%mitochondria

% spatial positions at which glucose is evaluated
% index 1 = point on domain edge
xpos = linspace(0,Lh,opt.gpts)'; 

if (~isempty(opt.startgluc))
    if (length(opt.startgluc) ~= opt.gpts); error('starting distrib has wrong size'); end
    gluc_init = opt.startgluc;
else
    lmdh = sqrt(Dh./(kgh*opt.nmito*Lh)); %lambda-hat
    gluc_init =  c0h * cosh((xpos-Lh/2)./(Lh * lmdh)) ./ cosh(0.5/lmdh);
end
gluc = gluc_init;
d2g = zeros(opt.gpts+2,1); %increased size to cater for new end points
dtg = zeros(opt.gpts+2,1); %increased size to cater to new end points
Tmito = zeros(opt.gpts+2,1);
gluc_new = zeros(opt.gpts+2,1);
Smito = zeros(opt.gpts+2,1);
ksx = zeros(opt.gpts+2,1);

ftc = 0; %flag for failing to converge. Is 1 when fails to converge. 
normdtg = inf;
dtcutoff = opt.dttol;%*(kgh*opt.Km/(opt.Km+1));
spacing = dx; %integration spacing

%% Iterative process
%continues till steady state
%steady state condition set by time derivative being small enough
step = 0;
while (normdtg > dtcutoff)
    %define glucose concentrations before x0 and after xL
    %to fix derivatives in terms of permeability at the edges
    gluc_start = gluc(2) + (2*Ph*(c0h - gluc(1))*dx / Dh);
    gluc_end = gluc(opt.gpts-1) + (2*Ph*(c0h - gluc(opt.gpts))*dx / Dh);
    
    %define gluc_new to be the new glucose vector to work with
    gluc_new = [gluc_start; gluc; gluc_end];    

    %Calculate distribution of total number of mitochondria
    ksx = ksh * opt.Km * gluc_new ./ (opt.Km + gluc_new);
    ksx_int = spacing * trapz(ksx);
    Tmito = (ksx/kwh + 1) ./ (Lh + (ksx_int/kwh));
    Smito = (ksx/kwh) ./ (Lh + (ksx_int/kwh));
    Smito_int = spacing * trapz(Smito);
    
    %Calculate the change in glucose concentration
    d2g(2:end-1) = (gluc_new(3:end)+gluc_new(1:end-2) - 2*gluc_new(2:end-1))/dx^2; %space double derivative
    % time derivative of glucose
    dtg(2:end-1) = Dh*d2g(2:end-1) - (kgh * opt.Km * opt.nmito) * (gluc_new(2:end-1) .* Tmito(2:end-1)) ./ (opt.Km + gluc_new(2:end-1));
    normdtg = norm(dtg);
    gluc_new = gluc_new+dtg*opt.delt;
       
    if (any(gluc_new < -1e-3))
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
        figure(3)
        plot(xpos,gluc_init,'k--')
        hold all
        plot(xpos,gluc,'b.-')        
        plot(xpos,Tmito(2:end-1)*1000,'r.-')
    title(sprintf('Step %d, normdtg= %f',step,normdtg))
        hold off
        drawnow
    end
    gluc = gluc_new(2:end-1);
end
end