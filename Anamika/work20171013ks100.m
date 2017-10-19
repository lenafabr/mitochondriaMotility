%% Run wrapper for high ks (ks = 100)
%changing lambda_hat through Kg
%changing c0_hat = c0/Km
 

% set up parameter values different from default

l_llim = -1.75;
l_ulim = 1;
nl = 101;
c0_llim = -3;
c0_ulim = 2;
nc0 = 102;
options.gpts = 100;
options.nmito = 500;
options.L = 500;
options.msize = 1;
options.D = 140;
options.dodisplay = 0;
options.nmito = 100;
options.dttol = 1e-3;
options.nsteps = 1e4;
options.ks = 100;


%run the function
%change lambda hat at every iteration, thereby changing kg
%change c0 at each iteration, thereby changing A

c0list = logspace(c0_llim,c0_ulim,nc0);
llist = logspace(l_llim,l_ulim,nl);

gluc_all = zeros(options.gpts,nl,nc0);
Tmito_all = zeros(options.gpts,nl,nc0);

for i = 1:1:nl
    lambda_hat(i) = llist(i);
    options.kg = options.D ./ (options.nmito * options.msize * options.L * (lambda_hat(i)^2));
    for j = 1:1:nc0
        options.c0 = c0list(j);
        options.cend = options.c0;
        [gluc,Tmito,Smito,Smito_int,normdtg,gluc_init,opt,xpos,lmdh,ftc] = runiterativesims(options);
        ftc_matrix(i,j) = ftc;
        option_list(i,j) = opt;
        gluc_init_all(:,i,j) = gluc_init;
        gluc_all(:,i,j) = gluc;
        Tmito_all(:,i,j) = Tmito;
        Smito_int_all(:,i,j) = Smito_int;
        Smito_all(:,i,j) = Smito;
        var_mito(i,j) = var(xpos,Tmito) ; %variance in mitochondria position distribution;
        varmetric(i,j) = 6*var_mito(i,j)/options.L^2 - 0.5;
    end
end
%save the workspace
formatOut = 'yyyymmdd';
date = datestr(datetime('today'),formatOut);
%save workspace with today's date'
filename = strcat('workspace_',date,'ks100');
save (filename);

%% get surface plot of varmetric vs c0 and lambda
figure;
varmetric = 6*var_mito/options.L^2 - 0.5;
colormap jet;
pcolor(log10(c0list),log10(llist),varmetric); shading flat
xlabel('log10(c0)')
ylabel('log10(lambda)')
title(sprintf('varmetric plot for high ks'));
%% look at conc necessary to achieve variance cutoff 
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

%%
% plot
% ncL = size(c0cutoffL);
% 
% c0cutoffLplot(ncL(2)+1:nl) = 0;
% c0cutoffLplot(1:ncL(2)) = c0cutoffL(1,:);
figure(2);
ll = logspace(-2,-1);
%loglog(lambda_hat,c0cutoffU,'r',lambda_hat,c0cutoffL, 'r','LineWidth',2)%,ll,1./ll.^3,ll,50./ll.^2,ll,30./ll)
loglog(lambda_hat,c0cutoffU,'r','LineWidth',2)%,ll,1./ll.^3,ll,50./ll.^2,ll,30./ll)
xlim([1e-2,1])
%
%loglog(lambda_hat,c0cutoff(:,A2_ind),lambda_hat,0.07./lambda_hat.^2)
xlabel('lambda hat')
ylabel(sprintf('conc to get %f cutoff',cutoff))
