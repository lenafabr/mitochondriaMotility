function [gluc, mitopos, mitostate, opt] = runmitosim(options)
% 
% %% options for testing
% options = struct();
% options.nmito = 10;
% 
% options.L = 100;
% options.D = 1e3;
% %options.startpos = 50;
% options.pstartwalk = 1;
% options.nstep = 1e5;
% options.restart = 1;
% options.D = 1e3;
% 
% options.showevery = 100;
% 
% options.cend=0.1;
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
opt.kg = 0.2; % rate of glucose consumption (if linear) 
opt.Km = 1;



% starting glucose distribution
% default is to start linear
opt.startgluc = [];

% fix permanent glucose distribution; do not evolve it
opt.fixgluc = [];

opt.nmito = 1; % number of mitochondria
opt.gpts = 100; % number of discrete spatial points for evaluating gluc concentration
opt.delt = 1e-4; % time-step
opt.nstep = 1000; % number of steps to run

% boundary conditions on the far size
% positive = fixed concentration at the boundary
% negative = reflecting boundary
opt.cend = 1; 

% starting position of mitochondria
% default (<0) means start uniformly
opt.startpos = -1;
% prob mitochondria start in walking state
% default is use equilibrium probability at c0 conc
opt.pstartwalk = opt.kw/(opt.kw + opt.ks*opt.c0);

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
% velh = opt.vel/opt.msize/opt.kg;
% kwh = opt.kw/opt.kg;
% ksh = opt.ks/opt.kg*opt.c0;
% Dh = opt.D/opt.msize^2/opt.kg;
velh = 1;
kwh = opt.kw/opt.vel*opt.msize
ksh = opt.ks/opt.vel*opt.c0*opt.msize;
Dh = opt.D/opt.msize/opt.vel;
kgh = opt.kg*opt.msize/opt.vel;
Kmh = opt.Km/opt.c0;

% spatial resolution
dx = Lh/(opt.gpts - 1);

%%
% initialize variables


% spatial positions at which glucose is evaluated
xpos = linspace(0,Lh,opt.gpts)';

if (opt.restart)
    % mitochondria center positions
    if opt.startpos>0
        % start at specific position
        mitopos = opt.startpos*ones(opt.nmito,1);
    else
        % start uniformly distributed btwn 1 and Lh-1
        mitopos = rand(opt.nmito,1)*(Lh-2)+1;
    end
    
    % glucose concentration
    % start with linear conc profile from c0 to cend
    if opt.cend<0
        error('reflecting boundary not yet implemented')
    else
        if (isempty(opt.startgluc))
            gluc = linspace(opt.c0,opt.cend,opt.gpts)';
        else
            if (length(opt.startgluc) ~= opt.gpts); error('starting distrib has wrong size'); end
            gluc = opt.startgluc;
        end
    end
    
    
    % state of the mitochondria (walking or not)
    % each mitochondria starts walking with probability p
    u = rand(opt.nmito,1);
    mitostate = u<=opt.pstartwalk;
    
    % set initial direction of walking (randomly)
    u = rand(opt.nmito,1);
    mitostate = ((u<=0.5)*2-1).*mitostate;
end

if (~isempty(opt.fixgluc))
    if (length(opt.fixgluc) ~= length(xpos)); error('wrong length of fixgluc'); end
    gluc = opt.fixgluc;
end

mitopos0 = mitopos;
%% evolve the system over time
d2g = zeros(opt.gpts,1);
dtg = d2g;

% probability of starting on each timestep
pstart = 1 - exp(-kwh*opt.delt);

if (opt.restart); curtime = 0; end

for step = 1:opt.nstep
    
    if (any(mitopos>Lh-0.5) | any(mitopos<0.5))
        error('bad mito positions')
    end
    
    if (isempty(opt.fixgluc))
    % ---------
    % evolve forward the glucose concentration by 1 time step
    % ----------
    % second derivative of gluc
    d2g(2:end-1) = (gluc(3:end)+gluc(1:end-2) - 2*gluc(2:end-1))/dx^2;
    % time derivative due to diffusion
    dtg = Dh*d2g;
    % time derivative due to glucose consumption
    for mc = 1:opt.nmito
        % get indices of spatial points within this mitochondria
        ind1 = floor((mitopos(mc)-0.5)/dx)+1;
        ind2 = floor((mitopos(mc)+0.5)/dx);        
        % dimensionless consumption rate of 1
        % note: overlapping mitochondria will consume twice as fast
        dtg(ind1:ind2) = dtg(ind1:ind2)-kgh*Kmh*gluc(ind1:ind2)./(Kmh+gluc(ind1:ind2));
    end
    
    % fix boundary conditions
    dtg(1) = 0; 
    if (opt.cend>0)
        % fixed conc at far end
        dtg(end) = 0;
    else
        % reflecting boundary at far end
        error('reflecting boundary not yet implemented')
    end
    
    % evolve forward
    gluc = gluc+dtg*opt.delt;    
    
    if (any(gluc<-1e-3))
        error('negative concentrations!')
    end    
    end
    
    % move the walking mitochondria
    walkind = find(mitostate);
    stopind = find(~mitostate);
    mitopos(walkind) = mitopos(walkind) + velh*opt.delt*mitostate(walkind);
    
    %% reflect mitochondria back if hitting the boundary
    for mc = walkind'
        if (mitopos(mc)<0.5)
            mitostate(mc)=1;
            mitopos(mc) = 0.5+(0.5-mitopos(mc));
        elseif mitopos(mc)>Lh-0.5
            mitostate(mc) = -1;
            mitopos(mc) = Lh-0.5 - (mitopos(mc)-Lh+0.5);
        end
    end
  
    %%
    % decide which mitochondria stop
    % glucose concentrations at mitochondria positions
    glucmito = interp1(xpos,gluc,mitopos(walkind));
    stoprate = ksh*Kmh*glucmito/(Kmh+glucmito);
    pstop = 1-exp(-stoprate*glucmito*opt.delt);
    u = rand(length(walkind),1);
    mitostate(walkind) = mitostate(walkind).*(1 - (u<=pstop));
       
    % decide which mitochondria start after this step   
    u = rand(length(stopind),1);
    restartind = find(u<=pstart);
    mitostate(stopind(restartind)) = (rand(length(restartind),1)<=0.5)*2 - 1; % start in random direction
    
    if (opt.dodisplay & mod(step,opt.showevery)==0)
        % plot glucose concentration
        plot(xpos,gluc,'.-')
        % plot mitochondria positions
        hold all
        ymin = min(gluc); ymax = max(gluc);
        for mc = 1:opt.nmito
            if (mitostate(mc)==0)
                %plot([mitopos(mc) mitopos(mc)], [ymin,ymax],'r','LineWidth',2)
                plot([mitopos(mc)], [0.5],'or','LineWidth',2)
            else
                %plot([mitopos(mc) mitopos(mc)], [ymin,ymax],'LineWidth',2,'Color',[0,0.5,0])
                 plot([mitopos(mc)], [0.5],'o','Color',[0,0.5,0],'LineWidth',2)
            end
            mval = (var(mitopos)-opt.L^2/12)/(opt.L^2/6);
            title(sprintf('Step %d: variance metric = %0.3f',step,mval))
            %set(gca,'FontSize',16)
            %legend('glucose','stopped mito', 'walking mito')
        end
        hold off
        %ylim([0,opt.c0*1.5])
        %legend('glucose','stopped mito', 'walking mito')
        drawnow
    end
    
    curtime = curtime + opt.delt;
end