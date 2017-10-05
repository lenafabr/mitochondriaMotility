% run iterative sim (check glucose profile and mitochondria distrib)

options.ks = 5; 
options.kg = 1;
options.c0 = 1;
options.cend = options.c0;
options.nmito = 100;
nitr = 10;
options.dodisplay=1;
options.showevery=100;
options.delt = 0.001;
options.L = 500;

options.gpts = 500;
options.nstep = 1e5;
options.dttol = 1e-3;
[gluc,Tmito,Smito,Smito_int,normdtg,gluc_init,opt,xpos,lmdh,ftc] = runiterativesims(options);

varmito = var(xpos,Tmito);
varmetric = varmito*6/options.L^2 - 0.5

%%
% sample initial mitochondria positions from Tmito

% tmito values at histogram bin edges
distrib = (Tmito(2:end) + Tmito(1:end-1))/2;
cTmito = cumsum(distrib)*(xpos(2)-xpos(1));
cTmito = cTmito/cTmito(end);
%xposdistrib = (xpos(2:end) + xpos(1:end-1))/2;
%xposdistrib = [0; xposdistrib];
xposdistrib = xpos;
cTmito = [0;cTmito];
plot(xposdistrib,cTmito)

% sample from the cumulative distribution
clear startpos
nsamp = 70;%options.nmito;
uvals = rand(1,nsamp);
startpos = interp1(cTmito,xposdistrib,uvals);


startpos = max(startpos,0.5);
startpos = min(startpos,options.L-0.5);

% compare histograms
[freq,bins] = hist(startpos,50)
plot(bins,freq)
db = bins(2)-bins(1);
plot(bins,freq/nsamp/db)
hold all
plot(xpos,Tmito)
hold off

var(startpos)*6/options.L^2 - 0.5

%% see sampled mitochondria
plot(startpos,zeros(size(startpos)),'o')
%% run discrete sims for comparison
options.nstep = 2e5;
options.kw = 1;
options.showmito = 1;
options.showevery = 500;

% TRY STARTING UNIFORM
options.startpos = -1;
% TRY STARTING AT EDGES
%options.startpos = [10,490];
% TRY STARTING AT CONTINUOUS SOLUTION
%options.startpos = startpos';

options.delt = 0.05;
nitr = 100;
options.pstartwalk = 1;
options.startgluc = gluc;
%
clear varmito gluc_dis mitopos_dis
%%
for j = 2:1:nitr
    
    % starting positions distributed according to continuous results
%     uvals = rand(1,nsamp);
%     startpos = interp1(cTmito,xposdistrib,uvals);
%     startpos = max(startpos,0.5);
%     startpos = min(startpos,options.L-0.5);
%     options.startpos = startpos';


    [gluc, mitopos, mitostate, opt] = runmitosim_michaelis2(options);
    varmito(j) = var(mitopos) ; %variance in mitochondria position distribution;
    gluc_dis(j,:) = gluc;
    mitopos_dis(j,:) = mitopos;
    
    [j varmito(j)]
    
    save('/home/ekoslover/proj/mitochondriaMotility/results/discretesims_500mito_unif_ks5_20170928.mat')
end