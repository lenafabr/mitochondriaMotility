%pull out peak values rowwise and the corresponding values of A
%plot the c0 vs lambda_hat
A2_i = 1;
A2_f = 10;
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

A2_ind = 1;
c0vals = A*options.kw/ks(A2_ind);
varmetric = 6*var_mito/options.L^2 - 0.5;
pcolor(log10(c0vals),log10(lambda_hat),varmetric(:,:,A2_ind)); shading flat
xlabel('log10(c0)')
ylabel('log10(lambda-hat)')

%% look at conc necessary to achieve variance cutoff (upper end)
A2_ind = 3;
cutoff = 0.03;
c0vals = A*options.kw/ks(A2_ind);
for i = 1:1:length(lambda_hat)
    [M,I] = max(varmetric(i,:,A2_ind)');
    if (M<cutoff || varmetric(i,end,A2_ind)>cutoff)
        c0cutoffU(i,A2_ind) = NaN;
    else
        c0cutoffU(i,A2_ind) = interp1(varmetric(i,I:end,A2_ind),c0vals(I:end),cutoff);
    end
end


%% look at conc necessary to achieve variance cutoff (lower end)

c0vals = A*options.kw/ks(A2_ind);
for i = 1:1:length(lambda_hat)
    [M,I] = max(varmetric(i,:,A2_ind)');
    if (M<cutoff || varmetric(i,end,A2_ind)>cutoff)
        c0cutoffL(i,A2_ind) = NaN;
    else
        c0cutoffL(i,A2_ind) = interp1(varmetric(i,1:I,A2_ind),c0vals(1:I),cutoff);
    end
end

%%
ll = logspace(-2,-1);
loglog(lambda_hat,c0cutoffU, 'r', lambda_hat,c0cutoffL,'b',ll,0.07./ll.^2,ll,0.11./ll)
%%
loglog(lambda_hat,c0cutoff(:,A2_ind),lambda_hat,0.07./lambda_hat.^2)
xlabel('lambda hat')
ylabel('conc to get 0.03 cutoff')