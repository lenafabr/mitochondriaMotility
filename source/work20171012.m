load('workspace_20171006ksc0.mat')

%% plot varmetric vs c0 and kshat
varmetric = 6*var_mito/options.L^2 - 0.5;
figure(1)
colormap jet;
pcolor(log10(c0list),log10(kslist),varmetric); shading flat
xlabel('log10(c0)')
ylabel('log10(ks_hat)')
title(sprintf('varmetric plot for lamhat=0.14'));

%% plot fraction stopped vs c0 and kshat
figure(2)
colormap jet;
pcolor(log10(c0list),log10(kslist),Smito_int_plot); shading flat
xlabel('log10(c0)')
ylabel('log10(ks_hat)')
title(sprintf('fraction stopped for lamhat=0.14'));

