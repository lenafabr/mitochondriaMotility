%script to plot variation measure
%first load workspace
c0 = 0.01:0.01:1.1;
var_c0 = var(mitopos_c0,1);
var_l = 500^2/12;
var_delta = 500^2/4;
figure 
plot(c0,(var_c0-var_l)/(var_delta-var_l));
title('Relative variation plot for 50 itr, 10000 time steps')
xlabel('c0') % x-axis label
ylabel('Var(mitopos) - VarU / Var\delta - VarU') % y-axis label
str = 'Km = 0.1/c0,Ks = 200*c0,Kw = 0.0100,Kg = 2';
dim = [.5 .5 .3 .3];
annotation('textbox',dim,'String',str,'FitBoxToText','on');