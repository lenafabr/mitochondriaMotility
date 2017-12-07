%Generating ks vs c0 surface plot
%% ks and c0 - from old data
load ('workspace_20171006ksc0.mat')
colormap jet;
pcolor(log10(c0list),log10(kslist),varmetric); shading flat
xlabel('log10(c0)')
ylabel('log10(k_s)')
title(sprintf('Variance Metric for varying k_s and c0,\lambda = %',lambda_hat));
figure
colormap jet;
pcolor(log10(c0list),log10(kslist),Smito_int_plot); shading flat
xlabel('log10(c0)')
ylabel('log10(k_s)')
title(sprintf('Fraction of mitochondria stopped for varying k_s and c0, \lambda = %',lambda_hat));
%% ks and c0 - from new data, different dimensionless units follows
%% Run wrapper for changing ks and c0
%lambda_hat fixed at 0.14
%changing c0_hat = c0/Km
%changing ks
%nmito = 100

% set up parameter values different from default
lambda_hat = 0.14;
options.D = 140;
c0_llim = -2;
c0_ulim = 4;
nc0 = 102;
ks_llim = -2;
ks_ulim = 3;
nks = 101;
options.gpts = 100;
options.nmito = 100;
options.L = 500;
options.msize = 1;
options.kg = options.D ./ (options.nmito * options.msize * options.L * (lambda_hat^2));
options.dodisplay = 0;
options.dttol = 1e-3;
options.delt = 1e-5/2;
options.nstep = 3*1e5;

%run the function
%change lambda hat at every iteration, thereby changing kg
%change c0 at each iteration, thereby changing A

c0list = logspace(c0_llim,c0_ulim,nc0);
kslist = logspace(ks_llim,ks_ulim,nks);

Gstat_all = zeros(options.gpts,nc0,nks);
Tmito_all = zeros(options.gpts,nc0,nks);


for i = 1:1:nks
    options.ks = kslist(i);
    for j = 1:1:nc0
        options.c0 = c0list(j);
        options.cend = options.c0;
        [gluc,Tmito,Smito,Smito_int,normdtg,gluc_init,opt,xpos,lmdh,ftc] = runiterativesims(options);
        normdtg_matrix(i,j) = normdtg;
        ftc_matric(i,j) = ftc;
        option_list(i,j) = opt;
        gluc_init_all(:,i,j) = gluc_init;
        gluc_all(:,i,j) = gluc;
        Tmito_all(:,i,j) = Tmito;
        Smito_int_all(:,i,j) = Smito_int;
        Smito_all(:,i,j) = Smito;
        var_mito(i,j) = var(xpos,Tmito) ; %variance in mitochondria position distribution;
        varmetric(i,j) = 6*var_mito(i,j)/options.L^2 - 0.5; 
    end
    percent_completed = (i/nks * 100)
end

%save the workspace
formatOut = 'yyyymmdd';
date = datestr(datetime('today'),formatOut);
%save workspace with today's date'
filename = strcat('workspace_',date,'ksc0_newdim');
save (filename);
%% get surface plot of varmetric vs c0 and ks
figure;
varmetric = 6*var_mito/options.L^2 - 0.5;
colormap jet;
pcolor(log10(c0list),log10(kslist),varmetric); shading flat
xlabel('log10(c0)')
ylabel('log10(ks)')
title(sprintf('varmetric plot for lmdh = 0.14'));


%% get surface plot of Smito_int vs c0 and ks
figure;
colormap jet;
pcolor(log10(c0list),log10(kslist),Smito_int_plot); shading flat
xlabel('log10(c0)')
ylabel('log10(ks)')
title(sprintf('fraction of mitochondria stopped for lmdh = 0.14'));


%% look at conc necessary to achieve variance cutoff 
% upper end

cutoff = 0.2;
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
loglog(lambda_hat,c0cutoffU,'r',lambda_hat,c0cutoffL, 'r','LineWidth',2)%,ll,1./ll.^3,ll,50./ll.^2,ll,30./ll)
xlim([1e-2,1])
%
%loglog(lambda_hat,c0cutoff(:,A2_ind),lambda_hat,0.07./lambda_hat.^2)
xlabel('lambda hat')
ylabel(sprintf('conc to get %f cutoff',cutoff))


