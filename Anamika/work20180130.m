%% To produce surface plots from existing data sets

%% get surface plot of varmetric vs c0 and lh
%load ks = 1e4 data
load('workspace_20180131lhc0_ks_1_10_100_1e4_Km0_1.mat');
figure;
varmetric(:,:,4) = 6*var_mito(:,:,4) - 0.5; %Lh = 1;
colormap jet;
pcolor(log10(c0list),log10(lhlist),varmetric(:,:,4)); shading flat
xlabel('log10(c0)')
ylabel('log10(\lambda)')
title(sprintf('Variance Metric Plot for k_s = 10000'));




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