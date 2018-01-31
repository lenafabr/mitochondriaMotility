%% set up default simulation parameters
opt = struct();

opt.kg = 1; % rate of glucose consumption
opt.c0 = 0.1; % fixed glucose concentration
opt.msize = 1; % mitochondria size

opt.L = 500; % domain size
opt.D = 140;% glucose diffusion coefficient
opt.kw = 1; % rate of starting a walk
opt.ks = 1; % rate of stopping is ks*[gluc]
opt.Km = 0.1;

% starting glucose distribution
% default is to start linear
opt.startgluc = [];

opt.nmito = 75; % number of mitochondria
opt.gpts = 100; % number of discrete spatial points for evaluating gluc concentration
opt.delt = 1e-5; % time-step
opt.nstep = 1e6; % number of steps to run

% starting position of mitochondria
% default (<0) means start uniformly
opt.startpos = -1;
% tolerance for "small time derivative"
opt.dttol = 1e-3;

% displaying plots
opt.dodisplay = 1;
opt.showevery = 1000;

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
ksh = opt.ks*cscale*tscale;
Dh = opt.D*tscale/lscale^2;
kgh = opt.kg*tscale;
c0h = opt.c0/cscale;
msizeh = opt.msize/lscale;
Kmh = opt.Km/cscale;

% spatial resolution
dx = Lh/(opt.gpts - 1);

%% Initialize start glucose concentration with analytical solution
%Analytical solution obtained by assuming uniform distribution of
%mitochondria
% spatial positions at which glucose is evaluated
% index 1 = point on domain edge
xpos = linspace(0,Lh,opt.gpts)';
lmdh = sqrt(Dh./(kgh*opt.nmito*msizeh*Lh)); %lambda-hat
if (~isempty(opt.startgluc))
    gluc_init=opt.startgluc;
else
    gluc_init = c0h * cosh((xpos-Lh/2)./(Lh * lmdh)) ./ cosh(0.5/lmdh);
end
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

dtcutoff = opt.dttol;

%% Solve the DE w/o using inbuilt function
% so that we can set our own tolerance values
step = 0;
while (normdtg > dtcutoff)
    
    %Calculate distribution of total number of mitochondria
    integrand =  Kmh * gluc ./ (Kmh + gluc);
    integ = dx/2 * (integrand(1) + integrand(end) + 2*sum(integrand(2:end-1)));
    Tmito = (Kmh * gluc ./ (Kmh + gluc))/ integ;

   %Calculate the change in glucose concentration
    d2g(2:end-1) = (gluc(3:end)+gluc(1:end-2) - 2*gluc(2:end-1))/dx^2; %space double derivative
    % time derivative of glucose
    dtg(2:end-1) = Dh*d2g(2:end-1) - ((opt.nmito * msizeh * kgh /integ)  * (((Kmh * gluc(2:end-1)) ./ (Kmh + gluc(2:end-1))).^2));
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
        var_mito= var(xpos,Tmito) ; %variance in mitochondria position distribution;
        varmetric = 6*var_mito/options.L^2 - 0.5;

        plot(xpos,gluc_init,'k--')
        hold all
        plot(xpos,gluc,'b.-')        
        plot(xpos,Tmito*Lh/2,'r.-')
        hold off
        title(sprintf('Step %d, dtg: %f, varmetric: %f',step,normdtg,varmetric))
        drawnow
    end
    Gstat = gluc;
    ksx_stat = (ksh * Kmh * Gstat) ./ (Kmh + Gstat);
    
end
