%% code similar to work20170809 - except for some corrections
%Runs corrected version of runmitosim_michaelis2 and cend = c0
%overall code runs discrete simulations for parameter values 
%which correspond to lambda = 10^-1.5 
%and give varmetric values of 0.2,0.3,0.5,0.7
%to figure out c0 values for given varmetric
load('workspace_08_09_1e5.mat');
lmdh_ind = 22;
A2_val = 20; A2_ind = 3;
%% find indexes of c0 for required varmetric value
varmetric_val = 0.2;
I1 = find(abs(varmetric(lmdh_ind,:,A2_ind) - varmetric_val)<0.01);
varmetric_val = 0.3;
I2 = find(abs(varmetric(lmdh_ind,:,A2_ind) - varmetric_val)<0.01);
varmetric_val = 0.5;
I3 = find(abs(varmetric(lmdh_ind,:,A2_ind) - varmetric_val)<0.01);
varmetric_val = 0.5;
I4 = find(abs(varmetric(lmdh_ind,:,A2_ind) - varmetric_val)<0.01);

%% run discrete simulations for the given parameters

c0_runlist = [c0vals(I1(1)),c0vals(I2(1)),c0vals(I3(1)),c0vals(I4(1))];
options.kg = options.D ./ (options.nmito * options.msize * options.L * (lambda_hat(lmdh_ind)^2));
options.ks = A2_val * options.kw ./ options.Km;
nitr = 100;
for j = 1:1:nitr
    for i = 1:1:4
        options.c0 = c0_runlist(i);
        options.cend = options.c0;
        [gluc, mitopos, mitostate, opt] = runmitosim_michaelis2(options);
        varmito_dis(j,i,:) = var(mitopos) ; %variance in mitochondria position distribution;
        gluc_dis(j,i,:) = gluc;
        mitopos_dis(j,i,:) = mitopos;
    end
    percent_completed = (j/nitr) * 100
end

%% run discrete simulations for c0s in different parts of phase diagram
 %cut off value for variance is 0.25
varmetric_val = 0.25;
lmdh_val = 0.08238;
A2_val = 20;

%calculate lower threshold index
Il = find(abs(c0vals-0.2597) < 0.005);
Iu = find(abs(c0vals-4.663) < 0.1);
%pull out c0 vals from different parts of the phase diagram
c0_pd_list = [c0vals(Il-2),c0vals(Il),c0vals(round((Il+Iu)/2)), c0vals(Iu), c0vals(Iu+2)];
options.kg = options.D ./ (options.nmito * options.msize * options.L * (lmdh^2));
options.ks = A2_val * options.kw ./ options.Km;
nitr = 100;
for j = 1:1:nitr
    for i = 1:1:5
        options.c0 = c0_pd_list(i);
        options.cend = options.c0;
        [gluc, mitopos, mitostate, opt] = runmitosim_michaelis2(options);
        varmito_dis(j,i,:) = var(mitopos) ; %variance in mitochondria position distribution;
        gluc_dis(j,i,:) = gluc;
        mitopos_dis(j,i,:) = mitopos;
    end
    percent_completed = j/nitr * 100
end
