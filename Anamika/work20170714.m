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
options.km = 100;
options.c0 = 1;

%options.ks = (1/4.8*1e-6)*c0*(10^-3*6e23/1000/1e12*4^2);
options.kw = 1 * 0.01;
options.delt=5e-2;
opt.msize = 1; % mitochondria size
options.dodisplay=0;
options.showevery=1;



%run the function
%change ks at every iteration, thereby changing A 
%Also change A* Lambda hat in every iteration, and thus calculate opt.kg
for i = 1:1:200
    options.ks = 0.01*i;
    A_var(i) = opt.ks * opt.c0/opt.kw;
    for j = 1:1:201
        Al(j) = 0.1*j;
        lmdh = Al ./ A_var(i);
        options.kg = options.D ./ (options.nmito * options.msize * options.L * lmdh^2);
        [gluc,Tmito,normdtg,gluc_init,opt,xpos,lmdh] = runiterativesims(options);
        var_metric(i,j) = var(xpos,Tmito) ; %variance in mitochondria position distribution;
    end
end

