% Regenerate phase space plot 
% Boundary conc vs lmdh with Kshat = 20, P = 0.1
% Compare with fixed end conc case

% For P = 0.1 use workspace_20171122p_01_l_c0.mat
% For Fixed end conc use workspace_20171127fixedc0.mat

%% First plot permeability case
load('workspace_20171128p_01_l_c0.mat')
 pcolor(log10(c0list),log10(lambda_hat),varmetric(:,:)); shading flat
 xlabel('log10(c0)')
 ylabel('log10(lambda-hat)')
 title(sprintf('P=%f',options.P))
%%
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
figure
for i = 1:1:length(lambda_hat)
    [M,I] = max(varmetric(i,:)');
    if (M<cutoff || varmetric(i,end)>cutoff)
        c0cutoffL(i) = NaN;
    else
        c0cutoffL(i) = interp1(varmetric(i,1:I),c0vals(1:I),cutoff);
    end
end


% plot
loglog(lambda_hat,c0cutoffU, 'b', lambda_hat,c0cutoffL,'b','LineWidth',2)
xlim([1e-2,1])

xlabel('lambda hat')
ylabel(sprintf('conc to get %f cutoff',cutoff))
hold on;


%% Then plot fixed c0 case
clear
load('workspace_20171130fixedc0.mat')
% figure
 pcolor(log10(c0list),log10(lambda_hat),varmetric(:,:)); shading flat
 xlabel('log10(c0)')
 ylabel('log10(lambda-hat)')
 title(sprintf('fixed c0'))
 %%
cutoff = 0.15;
%figure
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
loglog(lambda_hat,c0cutoffU, 'r', lambda_hat,c0cutoffL,'r','LineWidth',2)
xlim([1e-2,1])

xlabel('lambda hat')
ylabel(sprintf('conc to get %f cutoff',cutoff))
hold all;
