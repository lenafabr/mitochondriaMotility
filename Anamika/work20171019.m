load('workspace_20171012gstat2.mat')
%%
figure(1)
varmetric = 6*var_mito/options.L^2 - 0.5;
colormap jet;
pcolor(log10(c0list),log10(lambda_hat),varmetric); shading flat
xlabel('log10(c0)')
ylabel('log10(lambda)')
title(sprintf('varmetric plot for high ks'));

%%
load('workspace_08_09_1e6.mat')
varmetric=varmetric(:,:,1);
%% plot variance metric
%figure(2)
A2_ind =1;
ks = A2list(A2_ind)*options.kw/options.Km;
varmetric = 6*var_mito/options.L^2 - 0.5;
pcolor(log10(c0list),log10(lambda_hat),varmetric(:,:,A2_ind)); shading flat
xlabel('log10(c0)')
ylabel('log10(lambda-hat)')
title(sprintf('A2=%f',A2list(A2_ind)))

% --------------------
%% look at higher ks
load('workspace_20171013ks100.mat')

%% look at conc necessary to achieve variance cutoff (upper and lower)
% upper end
clear c0cutoffU c0cutoffL
cutoff = 0.15;
c0vals = logspace(c0_llim,c0_ulim,nc0);
for i = 1:1:length(lambda_hat)
    [M,I] = max(varmetric(i,:)');
    if (M<cutoff || varmetric(i,end)>cutoff)
        c0cutoffU(i) = NaN;
    else
        c0cutoffU(i) = interp1(varmetric(i,I:end),c0vals(I:end),cutoff);
    end
end


% (lower end)
c0vals = logspace(c0_llim,c0_ulim,nc0);
for i = 1:1:length(lambda_hat)
    [M,I] = max(varmetric(i,:)');
    if (M<cutoff || varmetric(i,end)>cutoff)
        c0cutoffL(i) = NaN;
    else
        c0cutoffL(i) = interp1(varmetric(i,1:I),c0vals(1:I),cutoff);
    end
end

% plot
ll = logspace(-2,-1);
loglog(lambda_hat,c0cutoffU, 'b', lambda_hat,c0cutoffL,'b','LineWidth',2)%,ll,1./ll.^3,ll,50./ll.^2,ll,30./ll)
xlim([1e-2,1])
xlabel('lambda hat')
ylabel(sprintf('conc to get %f cutoff',cutoff))

%% cutoff conc, upper end only
% upper end

cutoff = 0.15;
c0vals = logspace(c0_llim,c0_ulim,nc0);
for i = 1:1:length(lambda_hat)
    [M,I] = max(varmetric(i,:)');
    if (M<cutoff || varmetric(i,end)>cutoff)
        c0cutoffU(i) = NaN;
    else
        c0cutoffU(i) = interp1(varmetric(i,I:end),c0vals(I:end),cutoff);
    end
end
ll = logspace(-2,-1);
%loglog(lambda_hat,c0cutoffU,'r',lambda_hat,c0cutoffL, 'r','LineWidth',2)%,ll,1./ll.^3,ll,50./ll.^2,ll,30./ll)
loglog(lambda_hat,c0cutoffU,'r','LineWidth',2)%,ll,1./ll.^3,ll,50./ll.^2,ll,30./ll)
xlim([1e-2,1])
%
%loglog(lambda_hat,c0cutoff(:,A2_ind),lambda_hat,0.07./lambda_hat.^2)
xlabel('lambda hat')
ylabel(sprintf('conc to get %f cutoff',cutoff))

%% rerun high ks limit with 70 mitochondria
%% Run wrapper for run numeric sims for high ks
%changing lambda_hat through Kg
%changing c0_hat = c0/Km
%limit of high ks 

% set up parameter values different from default

l_llim = -1.75;
l_ulim = 1;
nl = 101;
c0_llim = -3;
c0_ulim = 2;
nc0 = 102;
options.gpts = 100;
options.L = 500;
options.msize = 1;
options.D = 140;
options.dodisplay = 1;
options.nmito = 70;
options.dttol = 1e-3;
options.nsteps = 1e4;


%run the function
%change lambda hat at every iteration, thereby changing kg
%change c0 at each iteration, thereby changing A

c0list = logspace(c0_llim,c0_ulim,nc0);
llist = logspace(l_llim,l_ulim,nl);

Gstat_all = zeros(options.gpts,nl,nc0);
Tmito_all = zeros(options.gpts,nl,nc0);

for i = 1:1:nl
    lambda_hat(i) = llist(i);
    options.kg = options.D ./ (options.nmito * options.msize * options.L * (lambda_hat(i)^2));
    for j = 1:1:nc0
        [i j nl nc0]
        options.c0 = c0list(j);
        options.cend = options.c0;
        [Gstat,Tmito,ksx_stat,gluc_init,opt,xpos,ftc] = gstatsim(options);
        ftc_matrix(i,j) = ftc;
        option_list(i,j) = opt;
        gluc_init_all(:,i,j) = gluc_init;
        Gstat_all(:,i,j) = Gstat;
        Tmito_all(:,i,j) = Tmito;
        var_mito(i,j) = var(xpos,Tmito) ; %variance in mitochondria position distribution;
        varmetric(i,j) = 6*var_mito(i,j)/options.L^2 - 0.5;
    end
end
%save the workspace
formatOut = 'yyyymmdd';
date = datestr(datetime('today'),formatOut);
%save workspace with today's date'
filename = strcat('workspace_',date,'gstat2_nmito70');
save (filename);

% -----------
%% Run iterative sims with high ks
%% iterative sims for same parameters
options.ks = 20; 
lh = 10^(-1.4);
kg = 140/(1*
options.kg = optdisc.kg;
options.c0 = optdisc.c0*10;
%options.ks = optdisc.ks; 
%options.kg = optdisc.kg;
%options.c0 = optdisc.c0;
options.cend = options.c0;
options.nmito = 500;%optdisc.nmito;

options.dodisplay=1;
options.showevery=100;
options.delt = 0.005;
options.L = 500;

options.gpts = 200;
options.nstep = 1e5;
options.dttol = 1e-3;
[gluc,Tmito,Smito,Smito_int,normdtg,gluc_init,opt,xpos,lmdh,ftc] = runiterativesims(options);

varmito = var(xpos,Tmito);
varmetric = varmito*6/options.L^2 - 0.5

% ---------------------
%% look at permeability=0, reflecting boundaries
options = struct();
options.kg = 0;
options.P = 0;
options.gpts=100;
options.delt = 1e-2; % time-step
options.nstep = 1e5; % number of steps to run
options.showevery=100;

lmdh=0.1;
Lh = 1;
xpos = linspace(0,Lh,options.gpts)'; 
gluc_init =  cosh((xpos-Lh/2)./(Lh * lmdh)) ./ cosh(0.5/lmdh);

plot(xpos,gluc_init,'.-')

dx=xpos(2)-xpos(1);
integ_init =dx*trapz(gluc_init)

options.startgluc = gluc_init;
%%
[gluc,Tmito,Smito,Smito_int,normdtg,gluc_init,opt,xpos] = permeablesims(options)
%
integ_final =dx*trapz(gluc)

% -----------
%% look at surface plots with finite permeability

% varying lambda-hat, c0-hat, with P=0.1
load('workspace_20171019perm_l_c0.mat')
%%
varmetric = 6*var_mito/options.L^2 - 0.5;
colormap jet;
pcolor(log10(c0list),log10(lambda_hat),varmetric); shading flat
caxis([0,0.5])
xlabel('log10(c0)')
ylabel('log10(lambda)')
title(sprintf('varmetric plot'));

%%
% varying lambda-hat, c0-hat, with P=0.1
load('workspace_20171018perm_ks_c0.mat')
%%
varmetric = 6*var_mito/options.L^2 - 0.5;
colormap jet;
pcolor(log10(c0list),log10(kslist),varmetric); shading flat
caxis([0,0.5])
xlabel('log10(c0)')
ylabel('log10(ks-hat)')
title(sprintf('varmetric plot'));

