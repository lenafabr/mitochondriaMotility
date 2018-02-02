%% To produce surface plots from existing data sets

%% get surface plot of varmetric vs c0 and lh
%load ks = {1,10,100,1e4} data
load('workspace_20180131lhc0_ks_1_10_100_1e4_Km0_1.mat'); %made using work20180116_2
figure;
varmetric = 6*var_mito- 0.5; %Lh = 1;
ksind = 4;
colormap jet;
pcolor(log10(c0list),log10(lhlist),varmetric(:,:,ksind)); shading flat
xlabel('log10(c0)')
ylabel('log10(\lambda)')
title(sprintf('Variance Metric Plot for k_s = %k',kslist(ksind)));

%% get surface plot from a different set of data
%load ks = {20,50,1000} data
%seems that ks = 1,10 are too small
load('workspace_20180201lhc0_ks_20_50_1000_Km0_1.mat')
varmetric = 6*var_mito- 0.5; %Lh = 1;
colormap jet;
ksind = 3;
pcolor(log10(c0list),log10(lhlist),varmetric(:,:,ksind)); shading flat
xlabel('log10(c0)')
ylabel('log10(\lambda)')
title(sprintf('Variance Metric Plot for k_s = %d',kslist(ksind)));

%% look at conc necessary to achieve variance cutoff 

for k = 1:1:nks
    % upper end
    cutoff = 0.2;
    c0vals = logspace(c0_llim,c0_ulim,nc0);
    for i = 1:1:length(lhlist)
        [M,I] = max(varmetric(i,:,k));
        if (M<cutoff || varmetric(i,end,k)>cutoff)
            c0cutoffU(i,k) = NaN;
        else
            c0cutoffU(i,k) = interp1(varmetric(i,I:end,k),c0vals(I:end),cutoff);
        end
    end
    
    
    % (lower end)
    c0vals = logspace(c0_llim,c0_ulim,nc0);
    for i = 1:1:length(lhlist)
        [M,I] = max(varmetric(i,:,k));
        if (M<cutoff || varmetric(i,end,k)>cutoff)
            c0cutoffL(i,k) = NaN;
        else
            c0cutoffL(i,k) = interp1(varmetric(i,1:I,k),c0vals(1:I),cutoff);
        end
    end
    cmap = [0,0,1;0.7,0,0.5;1,0,0];
    
    
    % plot
    %ll = logspace(-2,-1);
    loglog(lhlist,c0cutoffU(:,k),lhlist,c0cutoffL(:,k),'Color',cmap(k,:),'LineWidth',2)%,ll,1./ll.^3,ll,50./ll.^2,ll,30./ll)
    xlim([1e-2,1])
    hold on;
end
%loglog(lambda_hat,c0cutoff(:,A2_ind),lambda_hat,0.07./lambda_hat.^2)
xlabel('lambda hat')
ylabel(sprintf('conc to get %f cutoff',cutoff))