%% code similar to work20170809 - except for some corrections
%Similar to work20171822 except nmito = 700
%need to run iterative sims for A2 = 20 changing nmito = 700
%Runs corrected version of runmitosim_michaelis2 and cend = c0
%overall code runs discrete simulations for parameter values 
%which correspond to lambda = 10^-1.5 
%and give varmetric values of 0.2,0.3,0.5,0.7
%to figure out c0 values for given varmetric
load('workspace_20170824A2_20');
lmdh_ind = 22;
A2_val = 20; A2_ind = 1;
%% find indexes of c0 for required varmetric value
varmetric_val = 0.2;
I1 = find(abs(varmetric(lmdh_ind,:,A2_ind) - varmetric_val)<0.01);
varmetric_val = 0.3;
I2 = find(abs(varmetric(lmdh_ind,:,A2_ind) - varmetric_val)<0.01);
varmetric_val = 0.5;
I3 = find(abs(varmetric(lmdh_ind,:,A2_ind) - varmetric_val)<0.01);


%% run discrete simulations for the given parameters
c0vals = c0list;
c0_runlist = [c0vals(I1(1)),c0vals(I2(1)),c0vals(I3(1))];
options.kg = options.D ./ (options.nmito * options.msize * options.L * (lambda_hat(lmdh_ind)^2));
options.ks = A2_val * options.kw ./ options.Km;
nitr = 100;
for j = 1:1:nitr
    for i = 1:1:3
        options.c0 = c0_runlist(i);
        options.cend = options.c0;
        [gluc, mitopos, mitostate, opt] = runmitosim_michaelis2(options);
        varmito_dis(j,i,:) = var(mitopos) ; %variance in mitochondria position distribution;
        gluc_dis(j,i,:) = gluc;
        mitopos_dis(j,i,:) = mitopos;
    end
    percent_completed = (j/nitr) * 100
end
%save the workspace
formatOut = 'yyyymmdd';
date = datestr(datetime('today'),formatOut);
%save workspace with today's date'
filename = strcat('workspace_',date,'dis_1');
save (filename);

%% run discrete simulations for c0s in different parts of phase diagram
%5 parts, below, on between, on, above the curves
 %cut off value for variance is 0.20
 %concentration thresholds derived from visualize20170810 from phase
 %diagram, by eyeballing

varmetric_val = 0.20;
lmdh_val = 0.08722;
A2_val = 20;

%calculate lower threshold index
Il = find(abs(c0vals-0.1641) < 0.005);
Iu = find(abs(c0vals-5.477) < 0.1);
%pull out c0 vals from different parts of the phase diagram
c0_pd_list = [c0vals(Il-2),c0vals(Il),c0vals(round((Il+Iu)/2)), c0vals(Iu), c0vals(Iu+2)];
options.kg = options.D ./ (options.nmito * options.msize * options.L * (lmdh^2));
options.ks = A2_val * options.kw ./ options.Km;
options.nmito = 700;
nitr = 100;
for j = 1:1:nitr
    for i = 1:1:5
        options.c0 = c0_pd_list(i);
        options.delt = 5*1e-5;                                            
        options.cend = options.c0;
        [gluc, mitopos, mitostate, opt] = runmitosim_michaelis2(options);
        varmito_dis(j,i,:) = var(mitopos) ; %variance in mitochondria position distribution;
        gluc_dis(j,i,:) = gluc;
        mitopos_dis(j,i,:) = mitopos;
    end
    percent_completed = j/nitr * 100
end
%save the workspace
formatOut = 'yyyymmdd';
date = datestr(datetime('today'),formatOut);
%save workspace with today's date'
filename = strcat('workspace_',date,'dis_2');
save (filename);


%% Compare varmetric for discrete sims to give given value of varmetric 
load('workspace_08_22_dis_1');
load('workspace_08_09_1e5.mat');

varmito_dis_mean = mean(varmito_dis,1);

varmetric_dis_1= 6*varmito_dis_mean/options.L^2 - 0.5;

%% Compare varmetric for discrete sims in different part of the phase diagram
load('workspace_08_22_dis_2');
load('workspace_08_09_1e5.mat');

varmito_dis_mean = mean(varmito_dis,1);

varmetric_dis_2 = 6*varmito_dis_mean/options.L^2 - 0.5;

