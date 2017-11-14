%%

%For c0hat=100, lamhat = 0.02, kshat = 20, run sims with
%some diffrent permeabilities and compare to the fixed concentration case; 
%make sure it approaches the fixed conc as permeability gets high. 

options.gpts = 100;
options.nmito = 70;
options.L = 500;
options.msize = 1;
options.D = 140;
options.dodisplay = 0;
options.dttol = 1e-2;
options.delt = 1e-3;
options.nstep = 1e6;
options.c0 = 100;
options.cend = options.c0;
options.ks = 20; %kshat = 20 and kshat = ks*km/kw
lambda_hat = 0.02;
Plist = [0.1 0.2 0.5 1 10 20 50 100 1e4];
options.kg = options.D ./ (options.nmito * options.msize * options.L * (lambda_hat^2));
for i = 1:1:size(Plist,2)-1
    options.P = Plist(i);
    [gluc,Tmito,Smito,Smito_int,normdtg,gluc_init,opt,xpos,lmdh,ftc] = permeablesims(options);
    ftc_matrix(i) = ftc;
    option_list(i) = opt;
    gluc_init_all(:,i) = gluc_init;
    gluc_all(:,i) = gluc;
    Tmito_all(:,i) = Tmito;
    Smito_all(:,i) = Smito;
    Smito_int_all(i) = Smito_int;
    var_mito(i) = var(xpos,Tmito(2:end-1)) ; %variance in mitochondria position distribution;
    varmetric(i) = 6*var_mito(i)/options.L^2 - 0.5;
end
options.P = Plist(end);
options.delt = 1e-4;
[gluc,Tmito,Smito,Smito_int,normdtg,gluc_init,opt,xpos,lmdh,ftc] = permeablesims(options);
ftc_matrix(end) = ftc;
option_list(end) = opt;
gluc_init_all(:,end) = gluc_init;
gluc_all(:,end) = gluc;
Tmito_all(:,end) = Tmito;
Smito_all(:,end) = Smito;
Smito_int_all(end) = Smito_int;
var_mito(end) = var(xpos,Tmito(2:end-1)) ; %variance in mitochondria position distribution;
varmetric(end) = 6*var_mito(i)/options.L^2 - 0.5;

formatOut = 'yyyymmdd';
date = datestr(datetime('today'),formatOut);
%save workspace with today's date'
filename = strcat('workspace_',date,'Prange_l_0_02');
save (filename);

%% compare glucose profiles to profiles with end perm and fixed boundaries no P
%plot for different values of P
figure
load('workspace_20171113Prange_l_0_02.mat');

plot(xpos,gluc_all)
hold all;
clear;
load('workspace_08_09_1e6.mat');
lambda_hat_i = 13; %index associated with lambda_hat = 0.02
c0_i = 85; %index associated with c0 = 100;
plot(xpos,gluc_all(:,lambda_hat_i,c0_i,1),'r--');
plot(xpos,gluc_all(:,lambda_hat_i,1,1),'r--');
hold all;
% clear;
% 


