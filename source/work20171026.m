Lh=500;
Dh = options.D;
kgh = options.kg
xpos = linspace(0,Lh,opt.gpts)';
lmdh = sqrt(Dh./(kgh*opt.nmito*Lh)); %lambda-hat
gluc_init =  cosh((xpos-Lh/2)./(Lh * lmdh)) ./ cosh(0.5/lmdh);


%% check effect of nmito
options=struct();
options.gpts = 100;
options.L = 500;
options.msize = 1;
options.D = 140;
options.dodisplay = 0;
options.nmito = 70;
options.dttol = 1e-3;
options.nsteps = 1e4;
options.ks = 20;
options.kw = 1;
options.Km = 1;
options.P = 0.1;

lmdh =lambda_hat(2);
options.kg = options.D ./ (options.nmito * options.msize * options.L * (lmdh^2));
options.c0 = c0list(80);
options.cend = options.c0;

%%
[gluc,Tmito,Smito,Smito_int,normdtg,gluc_init,opt,xpos,lmdh,ftc] = permeablesims(options);
%
var_mito = var(xpos,Tmito(2:end-1))  %variance in mitochondria position distribution;
varmetric = 6*var_mito/options.L^2 - 0.5
%var_mito_70(2,80)
%% permeability sims
load('workspace_20171025p01_l_c0_70')
var_mito_70=var_mito;

load('workspace_20171019perm_l_c0.mat')
var_mito_100=var_mito;

%%
varmetric = 6*var_mito_100/options.L^2 - 0.5;
%%
pcolor(log10(c0list),log10(lambda_hat),varmetric); shading flat
xlabel('log10(c0)')
ylabel('log10(lambda-hat)')
colormap jet