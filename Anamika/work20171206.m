%Generating ks vs c0 surface plot
%% varmito vs ks and c0
colormap jet;
pcolor(log10(c0list),log10(kslist),varmetric); shading flat
xlabel('log10(c0)')
ylabel('log10(k_s)')
title(sprintf('Variance Metric for varying k_s and c0'));
%%