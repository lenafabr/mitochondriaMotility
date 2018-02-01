%% Run wrapper for changing c0 and lambda_hat
%ks fixed at 20,50,1000
%Km = 0.1
%changing c0_hat = c0/Km
%changing lambda_hat
%nmito = 75

% set up parameter values different from default
options.D = 140;
options.Km = 0.1;
c0_llim = -3.5;
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
options.dttol = 1e-2;
options.delt = 1e-5;
options.nstep = 1e7;

%run the function
%change lambda hat at every iteration, thereby changing kg
%change c0 at each iteration, thereby changing A

c0list = logspace(c0_llim,c0_ulim,nc0);
lhlist = logspace(lh_llim,lh_ulim,nlh);
kslist = [20,50,1000];
nks = 3;

Gstat_all = zeros(options.gpts,nc0,nlh,nks);
Tmito_all = zeros(options.gpts,nc0,nlh,nks);

for k = 1:1:nks
    options.ks = kslist(k);
    for i = 1:1:nlh
        lambda_hat = lhlist(i);
        options.kg = options.D ./ (options.nmito * options.msize * options.L * (lambda_hat^2));
        for j = 1:1:nc0
            options.c0 = c0list(j);
            options.cend = options.c0;
            [gluc,Tmito,Smito,Smito_int,normdtg,gluc_init,opt,xpos,lmdh,ftc] = runiterativesims(options);
            normdtg_matrix(i,j,k) = normdtg;
            ftc_matric(i,j,k) = ftc;
            option_list(i,j,k) = opt;
            gluc_init_all(:,i,j,k) = gluc_init;
            gluc_all(:,i,j,k) = gluc;
            Tmito_all(:,i,j,k) = Tmito;
            Smito_int_all(:,i,j,k) = Smito_int;
            Smito_all(:,i,j,k) = Smito;
            var_mito(i,j,k) = var(xpos,Tmito) ; %variance in mitochondria position distribution;
            varmetric(i,j,k) = 6*var_mito(i,j)- 0.5;
        end
        percent_completed = (i/nlh * 100) * (k/nks)
    end
end
%save the workspace
formatOut = 'yyyymmdd';
date = datestr(datetime('today'),formatOut);
%save workspace with today's date'
filename = strcat('workspace_',date,'lhc0_ks_20_50_1000_Km0_1');
save (filename);