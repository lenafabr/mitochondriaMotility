%run wrapper for runiterativesims - with log variation of c0,l (A,c0) & A2
%In a loop,lambda and kg are varied, and in the super loop A2 is varied
%A2 is varied by changing the ratio ks/kw. kw is kept fixed (at 1)
%define options structure

options = struct();
options.nmito = 14*5;

%c0 = 0.01

options.L = 500;
options.D = 140;
options.Km = 1;
%options.startpos = 50;
options.pstartwalk = 1;
options.nstep = 1e4;
options.restart = 1;
options.msize = 1;

options.kw = 1;
options.delt=5e-2;

options.dodisplay=0;
options.showevery=1;

ks_lim = 2;
nks = 100;
l_lim = 2;
nl = 101;
c0_lim = 2;
nc0 = 102;
%run the function
%change ks at every iteration, thereby changing A2
%change lambda hat at every iteration, thereby changing kg
%change c0 at each iteration, thereby changing A

for k = 1:1:nks
    logks = -ks_lim + (2*ks_lim/nks) * (k-1);
    options.ks = 10.^(logks);
    A2(k) = options.ks * options.Km ./ options.kw;
    for i = 1:1:nl
        log_lambda_hat = -l_lim + (2*l_lim/nl) * (i-1);
        lambda_hat(i) = 10 .^ (log_lambda_hat);
        options.kg = options.D ./ (options.nmito * options.msize * options.L * (lambda_hat(i)^2));
        for j = 1:1:nc0
            logc0 = -c0_lim + (2*c0_lim/nc0) * (j-1);
            options.c0 = 10 .^ (logc0);
            A(j) = options.ks * options.c0 / options.kw;
            [gluc,Tmito,normdtg,gluc_init,opt,xpos,lmdh] = runiterativesims(options)
            var_mito(i,j,k) = var(xpos,Tmito) ; %variance in mitochondria position distribution;
        end
    end
    percent_completed = k/nks * 100;
end
