%run wrapper for runiterativesims function
%define options structure

options = struct();
options.nmito = 14*5;

%c0 = 0.01

options.L = 500;
options.D = 140;
options.kg = 2;
%options.Km = 0.1/c0;
%options.startpos = 50;
options.pstartwalk = 1;
options.nstep = 1e4;
options.restart = 1;


%options.ks = (1/4.8*1e-6)*c0*(10^-3*6e23/1000/1e12*4^2);
options.kw = 1 * 0.01;
options.delt=5e-2;

options.dodisplay=1;
options.showevery=1;

%run the function
[gluc,Tmito,normdtg,gluc_init,opt] = runiterativesims(options)
plot(xpos,gluc_init,'k--')
hold all
plot(xpos,gluc,'b.-')
plot(xpos,Tmito*trapz(gluc_init),'r.-')
hold off
drawnow
