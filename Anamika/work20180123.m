%work file to generate plots

%% generate plots for fig 1(b) ie avg discrete distrib of mito

%first generate discrete histogram
load('workspace_20180119discretesims_100itr');
%concatenate mitopos_dis to one long array
mitodis = reshape(mitopos_dis,[1,size(mitopos_dis,1)*size(mitopos_dis,2)]);
hist_mito = histogram(mitodis,20,'Normalization','probability');
histogram(mitodis,100,'Normalization','probability');
hold on;

%then generate plot from iterative sims
options.nstep = 1e6;
options.Km = 0.1;
options.c0 = 0.1;
options.kw = 1;
options.L = 500;
options.ks=100;
options.kg=1;
options.dodisplay=0;
options.showevery=100;
options.nmito=75;

options.delt = 1e-5;
options.nstep = 1e7;
options.gpts = 100;
[gluc_itr,Tmito_itr,Smito_itr,Smito_int_itr,normdtg,gluc_init,opt,xpos,lmdh,ftc] = runiterativesims(options);

plot(xpos,Tmito_itr/100)

%% to plot figure 2 (Tmito for various c0) and the analytical soln
clear
% First plot Tmito from simulations
load('workspace_20180118c0_ks_200.mat')
T = [0,0,0.5;1,0,0];
x = [0;20];
cmap = interp1(x/20,T,linspace(0,1,20));
for c = 1:2:20
    plot(xpos,Tmito_all(:,c),'Color',cmap(c,:))
    hold all
end

c0 = options.Km;
options.kw = 1;
c0h = c0/options.Km;
ksh = options.ks * options.Km * (options.L)^2 / options.D;
Lh = 1;
gluc_calc = c0h * cosh((xpos-Lh/2)./(Lh * lmdh)) ./ cosh(0.5/lmdh);
dx = Lh/(options.gpts - 1);
ksx = ksh * options.Km * gluc_calc ./ (options.Km + gluc_calc);
ksx_int = dx * trapz(ksx);
kwh = options.kw * (options.L)^2 / options.D;
Tmito_calc = (ksx/kwh + 1) ./ (Lh + (ksx_int/kwh));
plot(xpos,Tmito_calc,'k--');

