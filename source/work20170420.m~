
%% options for testing
options = struct();
options.nmito = 14*5;

options.L = 500;
options.D = 200;
options.kg = 7e-0;
%options.startpos = 50;
options.pstartwalk = 1;
options.nstep = 1e5;
options.restart = 1;

options.showevery = 100;

options.cend=0.1;
options.ks = 1e-2;
options.kw = 1e-2;
options.delt=1e-2;

%options.startgluc = gluc;
%options.startgluc(end) = 0.1;

[gluc,mitopos,mitostate,opt] = runmitosim(options)
xpos = linspace(0,500,opt.gpts)';

%%
Km1d = 0.1e-3*6e23/1e3*(4e-4)^2*(1e-4)
kg = 2e5/Km1d