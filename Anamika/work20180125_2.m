%run iterative siumlations for several c0
%for ks = 100 for a new plot
%This is for the plot showing varying Tmito for c0, blue to red low c0 to
%high c0, with analytical solution shown 
%also want to get data set for glucose for a plot
lambda_hat = 0.06;
options.ks = 100;
options.D = 140;
c0_llim = -2;
c0_ulim = 2;
nc0 = 20;
options.D = 140;
options.msize = 1;
options.nmito = 75;
options.gpts = 100;
options.L = 500;
options.kg = 1;
options.Km = 0.1;

options.dodisplay = 0;
options.dttol = 1e-3;
options.delt = 1e-5;
options.nstep = 1e7;

%run the function
%change lambda hat at every iteration, thereby changing kg
%change c0 at each iteration, thereby changing A

c0list = logspace(c0_llim,c0_ulim,nc0);
gluc_all = zeros(options.gpts,nc0);
Tmito_all = zeros(options.gpts,nc0);



for j = 1:1:nc0
    options.c0 = c0list(j);
    options.cend = options.c0;
    [gluc,Tmito,Smito,Smito_int,normdtg,gluc_init,opt,xpos,lmdh,ftc] = runiterativesims(options);
    normdtg_matrix(j) = normdtg;
    ftc_matric(j) = ftc;
    option_list(j) = opt;
    gluc_init_all(:,j) = gluc_init;
    gluc_all(:,j) = gluc;
    gluc_over_c0(:,j) = gluc/(opt.c0/opt.Km); %dimensionless gluc over dimensionless c0
    Tmito_all(:,j) = Tmito;
    Smito_int_all(:,j) = Smito_int;
    Smito_all(:,j) = Smito;
    var_mito(j) = var(xpos,Tmito) ; %variance in mitochondria position distribution;
    varmetric(j) = 6*var_mito(j) - 0.5;
    percent_completed_c0 = (j/nc0 * 100)
end


%save the workspace
formatOut = 'yyyymmdd';
date = datestr(datetime('today'),formatOut);
%save workspace with today's date'
filename = strcat('workspace_',date,'c0_ks_100');
save (filename);
%% to plot figure 2 (Tmito for various c0) and the analytical soln
clear
% First plot Tmito from simulations
load('workspace_20180130c0_ks_100.mat')
T = [0,0,0.25;1,0,0];
x = [0;20];
cmap = interp1(x/20,T,linspace(0,1,20));
for c = 1:2:20
    plot(xpos,Tmito_all(:,c),'Color',cmap(c,:))
    hold all
end

c0 = 0.01; % 10 times lower than options.Km. Gives the extreme limit
options.kw = 1;
c0h = c0/options.Km;
ksh = options.ks * options.Km * (options.L)^2 / options.D;
Kmh = options.Km/options.Km;
Lh = 1;
gluc_calc = c0h * cosh((xpos-Lh/2)./(Lh * lmdh)) ./ cosh(0.5/lmdh);
dx = Lh/(options.gpts - 1);
ksx = ksh * Kmh * gluc_calc ./ (Kmh + gluc_calc);
ksx_int = dx * trapz(ksx);
kwh = options.kw * (options.L)^2 / options.D;
Tmito_calc = (ksx/kwh + 1) ./ (Lh + (ksx_int/kwh));
plot(xpos,Tmito_calc,'k--');

%plot gluc/c0 vs xpos 
figure
for c = 1:2:20
    plot(xpos,gluc_over_c0(:,c),'Color',cmap(c,:))
    hold all
end
plot(xpos,gluc_calc/c0h,'k--'); %analytical 
