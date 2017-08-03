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
options.nstep = 3*1e4;
options.restart = 1;
options.msize = 1;

options.kw = 1;
options.delt=5e-2;

options.dodisplay=0;
options.showevery=1;

A2_llim = -2;
A2_ulim = 2;
nA2 = 5;
l_llim = -2;
l_ulim = 0.5;
nl = 101;
c0_llim = -1;
c0_ulim = 2.5;
nc0 = 102;
%run the function
%change A2 at every iteration, thereby changing ks
%change lambda hat at every iteration, thereby changing kg
%change c0 at each iteration, thereby changing A

for k = 1:1:nA2
    logA2 = A2_llim + ((A2_ulim - A2_llim)/nA2) * (k-1);
    A2(k) = 10.^(logA2);
    options.ks = A2(k) * options.kw ./ options.Km;
    for i = 1:1:nl
        log_lambda_hat = l_llim + ((l_ulim - l_llim)/nl) * (i-1);
        lambda_hat(i) = 10 .^ (log_lambda_hat);
        options.kg = options.D ./ (options.nmito * options.msize * options.L * (lambda_hat(i)^2));
        for j = 1:1:nc0
            logc0 = c0_llim + ((c0_ulim - c0_llim)/nc0) * (j-1);
            options.c0 = 10 .^ (logc0);
            A(j) = options.ks * options.c0 / options.kw;
            [gluc,Tmito,Smito,Smito_int,normdtg,gluc_init,opt,xpos,lmdh,ftc] = runiterativesims(options);
            ftc_matrix(i,j,k) = ftc;
            gluc(:,i,j,k) = gluc;
            var_mito(i,j,k) = var(xpos,Tmito) ; %variance in mitochondria position distribution;
        end
    end
    percent_completed = (k/nA2)* 100
end
