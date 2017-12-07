lambda_hat = 0.2122;
c0 = 40.179;

%% Run with edge permeabilities, for a particular set of parameters
options.P = 0.1;
options.kg = options.D ./ (options.nmito * options.msize * options.L * (lambda_hat^2));
options.c0 = c0;
options.cend = options.c0;
options.dodisplay=1;
options.showevery=1000;
[gluc,Tmito,Smito,Smito_int,normdtg,gluc_init,opt,xpos,lmdh,ftc] = permeablesims(options);

var_mito = var(xpos,Tmito(2:end-1)) ; %variance in mitochondria position distribution;
varmetric = 6*var_mito/options.L^2 - 0.5;
    
%% run with fixed end concentrations
options.kg = options.D ./ (options.nmito * options.msize * options.L * (lambda_hat^2));
options.c0 = c0;
options.cend = options.c0;
options.dodisplay=1;
options.showevery=1000;
[gluc,Tmito,Smito,Smito_int,normdtg,gluc_init,opt,xpos,lmdh,ftc] = runiterativesims(options);

var_mito = var(xpos,Tmito) ; %variance in mitochondria position distribution;
varmetric = 6*var_mito/options.L^2 - 0.5

%% for debugging, run constant maxed-out consumption case
options.P = 0.1;
options.kg = options.D ./ (options.nmito * options.msize * options.L * (lambda_hat^2));
options.c0 = c0;
options.cend = options.c0;
options.dodisplay=1;
options.showevery=10000;
options.nstep=1e7;
options.gpts = 100;
[gluc,Tmito,Smito,Smito_int,normdtg,gluc_init,opt,xpos,lmdh,ftc] = permeablesims_debug(options);


%%
vel=1;
Lh = options.L/options.msize;
Dh = options.D/options.msize^2*options.msize/vel;
c0h = options.c0/options.Km;
Ph = options.P/options.msize*options.msize/vel;
kgh=options.kg*options.msize/vel;

A = kgh*options.nmito/(Lh*2*Dh);
B = c0h - A*Lh^2/4 - A*Lh*Dh/Ph;
solan = A*(xpos-options.L/2).^2 + B;
hold all
plot(xpos,solan,'-.')
hold off

dGdx1 = (solan(end)-solan(end-1))/(xpos(end)-xpos(end-1))
checkend = Dh*dGdx1 - Ph*(c0h-solan(end))
%
dGdx2 = (gluc(end)-gluc(end-1))/(xpos(end)-xpos(end-1));
checkend = Dh*dGdx2 - Ph*(c0h-gluc(end))

%%