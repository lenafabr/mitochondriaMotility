%% Run wrapper for changing c0 and lambda_hat
%ks fixed at 1, 10, 100
%get 1e4 data from workspace_lhc0newdim
%changing c0_hat = c0/Km
%changing lambda_hat
%nmito = 75

% set up parameter values different from default
lambda_hat = 0.06;
options.D = 140;
c0_llim = -2;
c0_ulim = 2;
nc0 = 102;
lh_llim = -2;
lh_ulim = 2;
nlh = 101;
options.D = 140;
options.msize = 1;
options.nmito = 75;
options.gpts = 100;
options.L = 500;

options.dodisplay = 0;
options.dttol = 1e-3;
options.delt = 1e-5;
options.nstep = 1e6;

%run the function
%change lambda hat at every iteration, thereby changing kg
%change c0 at each iteration, thereby changing A

c0list = logspace(c0_llim,c0_ulim,nc0);
lhlist = logspace(lh_llim,lh_ulim,nlh);
kslist = {1,10,100};

Gstat_all = zeros(options.gpts,nc0,nlh,3);
Tmito_all = zeros(options.gpts,nc0,nlh,3);

for k = 1:1:3
    options.ks = kslist{k};
    for i = 1:1:nlh
        lambda_hat = lhlist(i);
        options.kg = options.D ./ (options.nmito * options.msize * options.L * (lambda_hat^2));
        for j = 1:1:nc0
            options.c0 = c0list(j);
            options.cend = options.c0;
            [gluc,Tmito,Smito,Smito_int,normdtg,gluc_init,opt,xpos,lmdh,ftc] = runiterativesims(options);
            normdtg_matrix(i,j) = normdtg;
            ftc_matric(i,j) = ftc;
            option_list(i,j) = opt;
            gluc_init_all(:,i,j) = gluc_init;
            gluc_all(:,i,j) = gluc;
            Tmito_all(:,i,j) = Tmito;
            Smito_int_all(:,i,j) = Smito_int;
            Smito_all(:,i,j) = Smito;
            var_mito(i,j) = var(xpos,Tmito) ; %variance in mitochondria position distribution;
            varmetric(i,j) = 6*var_mito(i,j)/options.L^2 - 0.5;
        end
        percent_completed = (i/nlh * 100)
    end
end
%save the workspace
formatOut = 'yyyymmdd';
date = datestr(datetime('today'),formatOut);
%save workspace with today's date'
filename = strcat('workspace_',date,'lhc0_newdim_ks_1_10_100');
save (filename);