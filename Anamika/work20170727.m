%pull out peak values rowwise and the corresponding values of A
%plot the c0 vs lambda_hat
A2_i = 1;
A2_f = 36;
ks = A2*options.kw/options.Km;
for A2_ind = A2_i:1:A2_f
    [M,I] = max(var_mito(:,:,A2_ind)');
    for i = 1:1:101
        A_ind = I(i);
        c0_opt(i,A2_ind) = A(A_ind) * options.kw/ks(A2_ind);
    end
    %plot c0_opt
    plot(lambda_hat,c0_opt(:,A2_ind));
    hold on;
    
end


