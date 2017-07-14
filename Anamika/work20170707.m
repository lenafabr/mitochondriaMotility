%run wrapper for runiterativesims function - with changing parameters
%define options structure

options = struct();
options.nmito = 14*5;

%c0 = 0.01

options.L = 500;
options.D = 140;
%options.Km = 0.1/c0;
%options.startpos = 50;
options.pstartwalk = 1;
options.nstep = 1e4;
options.restart = 1;


%options.ks = (1/4.8*1e-6)*c0*(10^-3*6e23/1000/1e12*4^2);
options.kw = 1 * 0.01;
options.delt=5e-2;

options.dodisplay=0;
options.showevery=1;


%run the function
%change c0 at each iteration, thereby changing A
%change Kg at every iteration, thereby changing lambda hat
for i = 1:1:100
    options.kg = 0.1*i;
    for j = 1:1:100
        options.c0 = 0.1*j;
        [gluc,Tmito,normdtg,gluc_init,opt,lmdh] = runiterativesims(options)
        A_var(i,j) = opt.ks * opt.c0/opt.kw;
        A = A_var(i,j);
        Lambda_hat(i,j) = lmdh;
        Lh = opt.L/opt.msize;
        xpos = linspace(0,Lh,opt.gpts)';
        var_metric(i,j) = var(xpos,Tmito) ; %variance in mitochondria position distribution;
       
    end
end


% plot(xpos,gluc_init,'k--')
% hold all
% plot(xpos,gluc,'b.-')
% plot(xpos,Tmito*trapz(gluc_init),'r.-')
% hold off
% drawnow