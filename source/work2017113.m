options = struct();

options.L = 500;
options.Km = 1;
options.c0=50;
options.cend = 0;
options.kg = 100;
options.P = 0.1;
options.ks = 1000;
options.showevery=10;

options.gpts=200;
options.delt = 5e-3;

[gluc,Tmito,Smito,Smito_int,normdtg,gluc_init,opt,xpos,lmdh,ftc] = ...
    constantpsims(options);

%%
cext = linspace(options.c0,options.cend,length(xpos));
plot(xpos,options.P/(opt.nmito*options.kg/options.L)*cext,'g--')
hold all
plot(xpos,gluc,'b-')
hold off

%%
% analytical solution for uniform mito, linear kinetics
xvals = xpos-options.L/2;
keff = (options.kg*opt.nmito/options.L)+options.P;
beta = sqrt((keff+options.P)/opt.D);

lincoeff = (opt.cend-opt.c0)/opt.L;
constterm = (opt.cend+opt.c0)/2;

Gan = options.P/keff*(lincoeff*xvals+constterm) ...
    - options.P*lincoeff*sinh(beta*xvals)/(beta*keff*cosh(beta*options.L/2));

%plot(xpos,Gan,xpos, options.P/keff*(lincoeff*xvals+constterm),'--')

plot(xpos,gluc,'.-',xpos,Gan,'--')