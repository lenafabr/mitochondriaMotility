%% Ends permeable, run wrapper for changing c0 and lambda_hat and P
%ks fixed at 1e4
%changing c0_hat = c0/Km
%changing lambda_hat
%nmito = 75
%different permeability

% set up parameter values different from default
options.D = 140;
c0_llim = -2;
c0_ulim = 2;
nc0 = 102;
lh_llim = -2;
lh_ulim = 2;
nlh = 101;
options.Km = 0.1;
options.D = 140;
options.msize = 1;
options.nmito = 75;
options.gpts = 100;
options.L = 500;
options.ks = 200;

options.dodisplay = 0;
options.dttol = 1e-2;
options.delt = 5*1e-5;
options.nstep = 1e8;

%run the function
%change lambda hat at every iteration, thereby changing kg
%change c0 at each iteration, thereby changing A

c0list = logspace(c0_llim,c0_ulim,nc0);
lhlist = logspace(lh_llim,lh_ulim,nlh);
Plist = [0.01,0.05,0.1,1,10,100,1e4]; %0.05 is close to actual P values
nP = 7;

gluc_all = zeros(options.gpts,nlh,nc0,nP);
Tmito_all = zeros(options.gpts+2,nlh,nc0,nP);

for k = 1:1:7
    options.P = Plist(k);
    for i = 1:1:nlh
        lambda_hat = lhlist(i);
        options.kg = options.D ./ (options.nmito * options.msize * options.L * (lambda_hat^2));
        for j = 1:1:nc0
            options.c0 = c0list(j);
            options.cend = options.c0;
            [gluc,Tmito,Smito,Smito_int,normdtg,gluc_init,opt,xpos,lmdh,ftc] = permeablesims(options);
            normdtg_matrix(i,j,k) = normdtg;
            ftc_matric(i,j,k) = ftc;
            option_list(i,j,k) = opt;
            gluc_init_all(:,i,j,k) = gluc_init;
            gluc_all(:,i,j,k) = gluc;
            Tmito_all(:,i,j,k) = Tmito;
            Smito_int_all(:,i,j,k) = Smito_int;
            Smito_all(:,i,j,k) = Smito;
            var_mito(i,j,k) = var(xpos,Tmito(2:end-1)) ; %variance in mitochondria position distribution;
            varmetric(i,j,k) = 6*var_mito(i,j,k) - 0.5; %Lh = 1
        end
        percent_completed = (i/nlh * 100)*(k/7)
    end
end

%save the workspace
formatOut = 'yyyymmdd';
date = datestr(datetime('today'),formatOut);
%save workspace with today's date'
filename = strcat('workspace_',date,'endP_lh_c0');
save (filename);