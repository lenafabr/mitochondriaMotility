%pull out peak values rowwise and the corresponding values of A
%plot the c0 vs lambda_hat
A2_i = 1;
A2_f = 36;
ks = A2*options.kw/options.Km;
cmat = jet(A2_f-A2_i+1);
for A2_ind = A2_i:1:A2_f
    [M,I] = max(var_mito(:,:,A2_ind)');
    for i = 1:1:101
        A_ind = I(i);
        c0_opt(i,A2_ind) = A(A_ind) * options.kw/ks(A2_ind);
    end
    %plot c0_opt
    plot(log10(lambda_hat),log10(c0_opt(:,A2_ind)),'Color',cmat(A2_ind,:));
    hold on;
    
end

hold off

%% look at individual lam-hat vs c0 profiles
ks = A2*options.kw/options.Km;
A2_ind = 4;

c0vals = logspace(c0_llim,c0_ulim,nc0)

varmetric = 6*var_mito/options.L^2 - 0.5;
pcolor(log10(c0vals),log10(lambda_hat),varmetric(:,:,A2_ind)); shading flat
xlabel('log10(c0)')
ylabel('log10(lambda-hat)')
title(sprintf('A2=%f',A2(A2_ind)))

%% look at conc necessary to achieve variance cutoff (upper end)
cutoff = 0.22;
c0vals = logspace(c0_llim,c0_ulim,nc0);
for i = 1:1:length(lambda_hat)
    [M,I] = max(varmetric(i,:,A2_ind)');
    if (M<cutoff || varmetric(i,end,A2_ind)>cutoff)
        c0cutoffU(i) = NaN;
    else
        c0cutoffU(i) = interp1(varmetric(i,I:end,A2_ind),c0vals(I:end),cutoff);
    end
end


%% look at conc necessary to achieve variance cutoff (lower end)
c0vals = logspace(c0_llim,c0_ulim,nc0);
for i = 1:1:length(lambda_hat)
    [M,I] = max(varmetric(i,:,A2_ind)');
    if (M<cutoff || varmetric(i,end,A2_ind)>cutoff)
        c0cutoffL(i) = NaN;
    else
        c0cutoffL(i) = interp1(varmetric(i,1:I,A2_ind),c0vals(1:I),cutoff);
    end
end


%%
ll = logspace(-2,-1);
loglog(lambda_hat,c0cutoffU, 'c', lambda_hat,c0cutoffL,'c')%,ll,1./ll.^3,ll,50./ll.^2,ll,30./ll)
xlim([1e-2,1])
%
%loglog(lambda_hat,c0cutoff(:,A2_ind),lambda_hat,0.07./lambda_hat.^2)
xlabel('lambda hat')
ylabel('conc to get 0.25 cutoff')

%% look at glucose profiles
lind = 10;

c0want = 160;

c0vals = logspace(c0_llim,c0_ulim,nc0)

for A2_ind = [4,5]
    
    c0ind = round(interp1(c0vals,1:length(A),c0want));
    
    disp([A2_ind c0ind c0vals(c0ind)])
    
    plot(gluc(:,lind,c0ind,A2_ind))
    hold all
end
hold off


%% Rerun some sims to check gluc profiles

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

gluc_all = zeros(100,nA2);
Smito_all = zeros(100,nA2);
Smito_int_all = zeros(1,nA2);
clear var_mito
for k = [4,5]
    logA2 = A2_llim + ((A2_ulim - A2_llim)/nA2) * (k-1);
    A2(k) = 10.^(logA2);
    options.ks = A2(k) * options.kw ./ options.Km;
    
    i=10
        log_lambda_hat = l_llim + ((l_ulim - l_llim)/nl) * (i-1);
        lambda_hat(i) = 10 .^ (log_lambda_hat);
        options.kg = options.D ./ (options.nmito * options.msize * options.L * (lambda_hat(i)^2));
        
    
        options.c0 = 290;
        
        [gluc,Tmito,Smito,Smito_int,normdtg,gluc_init,opt,xpos,lmdh,ftc] = runiterativesims(options);
        
        ftc_matrix(k) = ftc;
        Smito_all(:,k) = Smito;
        Smito_int_all(k) = Smito_int;
        gluc_all(:,k) = gluc;
        var_mito(k) = var(xpos,Tmito) ; %variance in mitochondria position distribution;
     
end

%%
varmetric = 6*var_mito/options.L^2 - 0.5