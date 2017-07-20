%run wrapper for runiterativesims function - with changing parameters
%change Al and l logarithmically

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
options.km = 20;
options.c0 = 1;

%options.ks = (1/4.8*1e-6)*c0*(10^-3*6e23/1000/1e12*4^2);
options.kw = 1 * 0.01;
options.delt=5e-2;
options.msize = 1; % mitochondria size
options.dodisplay=0;
options.showevery=1;



%run the function
%change lambda hat in every iteration, calculate options.kg 
%Also change A* Lambda hat in every iteration, and thus calculate
%options.ks
nl = 200;
nAl = 201;
l_lim = 3; %lambda's absolute limit
Al_lim = 2; %Al's absolute limit
for i = 1:1:nl
    log_lambda_hat = l_lim - (2*l_lim/nl) * (i-1);
    lambda_hat(i) = 10 .^ (log_lambda_hat);
    for j = 1:1:201
        log_Al = -Al_lim + (2*Al_lim/nAl) * j;
        Al(j) = 10 .^ (log_Al);
        A = Al(j) ./ lambda_hat(i);
        options.kg = options.D ./ (options.nmito * options.msize * options.L * (lambda_hat(i)^2));
        options.ks = A * options.kw ./ options.c0;
        [gluc,Tmito,normdtg,gluc_init,opt,xpos,lmdh] = runiterativesims(options);
        var_metric(i,j) = var(xpos,Tmito) ; %variance in mitochondria position distribution;
    end
end