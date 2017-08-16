load('workspace_08_09_1e6.mat')

%%
load('workspace_08_09_1e5.mat')
%%
A2_ind =2;

ks = A2list(A2_ind)*options.kw/options.Km;
varmetric = 6*var_mito/options.L^2 - 0.5;
pcolor(log10(c0list),log10(lambda_hat),varmetric(:,:,A2_ind)); shading flat
xlabel('log10(c0)')
ylabel('log10(lambda-hat)')
title(sprintf('A2=%f',A2list(A2_ind)))

%% look at conc necessary to achieve variance cutoff 
% upper end
clear c0cutoffU c0cutoffL
cutoff = 0.1;
c0vals = logspace(c0_llim,c0_ulim,nc0);
for i = 1:1:length(lambda_hat)
    [M,I] = max(varmetric(i,:,A2_ind)');
    if (M<cutoff || varmetric(i,end,A2_ind)>cutoff)
        c0cutoffU(i) = NaN;
    else
        c0cutoffU(i) = interp1(varmetric(i,I:end,A2_ind),c0vals(I:end),cutoff);
    end
end


% (lower end)
c0vals = logspace(c0_llim,c0_ulim,nc0);
for i = 1:1:length(lambda_hat)
    [M,I] = max(varmetric(i,:,A2_ind)');
    if (M<cutoff || varmetric(i,end,A2_ind)>cutoff)
        c0cutoffL(i) = NaN;
    else
        c0cutoffL(i) = interp1(varmetric(i,1:I,A2_ind),c0vals(1:I),cutoff);
    end
end


% plot
ll = logspace(-2,-1);
loglog(lambda_hat,c0cutoffU, 'b', lambda_hat,c0cutoffL,'b')%,ll,1./ll.^3,ll,50./ll.^2,ll,30./ll)
xlim([1e-2,1])
%
%loglog(lambda_hat,c0cutoff(:,A2_ind),lambda_hat,0.07./lambda_hat.^2)
xlabel('lambda hat')
ylabel(sprintf('conc to get %f cutoff',cutoff))


% -------------
%% discrete sims for different variance metrics (what variance metric is reasonable to use?)
% generated with Anamika/work20170809.m
% discrete sims for varmetric values: 0.2, 0.3,0.5,0.7
load('workspace_08_09_dis_1.mat')
A2 = 20
lam = 10^-1.5
c0_runlist = [0.2605    0.5786    3.3477    3.3477];
mitopos_dis = mitopos_dis(1:4,:);


%% plot individual mitochondria
cc = 3;
subplot(1,2,1)
plot(mitopos_dis(cc,:),zeros(1,70),'o')
subplot(1,2,2)
hist(mitopos_dis(cc,:),20)

varmetric = 6*var(mitopos_dis(cc,:))/500^2 - 0.5

%%
load('workspace_08_09_dis_2.mat')