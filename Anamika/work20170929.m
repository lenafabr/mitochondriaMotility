%exact same as work20170829 - except for new set of parameters
%these new parameters give a more extreme varmetric value
%expect varmetric = .24 from iterative sims
% run iterative sim (check glucose profile and mitochondria distrib)

options.ks = 50; 
options.kg = 0.1;
options.c0 = 0.1;
options.cend = options.c0;
options.nmito = 500;
nitr = 10;
options.dodisplay=1;
options.showevery=100;
options.delt = 0.05;
options.L = 500;

[gluc,Tmito,Smito,Smito_int,normdtg,gluc_init,opt,xpos,lmdh,ftc] = runiterativesims(options);

varmito = var(xpos,Tmito);
varmetric = varmito*6/options.L^2 - 0.5

%%
% sample initial mitochondria positions from Tmito

% tmito values at histogram bin edges
distrib = (Tmito(2:end) + Tmito(1:end-1))/2;
cTmito = cumsum(distrib)*(xpos(2)-xpos(1));
cTmito = cTmito/cTmito(end);
xposdistrib = (xpos(2:end) + xpos(1:end-1))/2;
xposdistrib = [0; xposdistrib];
cTmito = [0;cTmito];
plot(xposdistrib,cTmito)

% sample from the cumulative distribution
clear startpos
nsamp = options.nmito;
for mc = 1:nsamp
    u = rand();
    startpos(mc) = interp1(cTmito,xposdistrib,u);
end

startpos = max(startpos,0.5);
startpos = min(startpos,options.L-0.5);

% compare histograms
[freq,bins] = hist(startpos,10)
plot(bins,freq)
db = bins(2)-bins(1);
plot(bins,freq/nsamp/db)
hold all
plot(xpos,Tmito)
hold off


%% run discrete sims for comparison
options.nstep = 1e6;
options.kw = 1;
options.showmito = 1;
options.showevery = 100;

% TRY STARTING UNIFORM
%options.startpos = -1;
% TRY STARTING AT EDGES
%options.startpos = [10,490];
% TRY STARTING AT CONTINUOUS SOLUTION
options.startpos = startpos';

options.delt = 0.05;
nitr = 100;
options.pstartwalk = 1;
options.startgluc = gluc;
%
clear varmito gluc_dis mitopos_dis
for j = 1:1:nitr
    [gluc, mitopos, mitostate, opt] = runmitosim_michaelis2(options);
    varmito(j) = var(mitopos) ; %variance in mitochondria position distribution;
    gluc_dis(j,:) = gluc;
    mitopos_dis(j,:) = mitopos;
    
    [j varmito(j)]
end
%save the workspace
formatOut = 'yyyymmdd';
date = datestr(datetime('today'),formatOut);
%save workspace with today's date'
filename = strcat('workspace_',date,'continuous_sol');
save (filename);
