%% figure out iterative vs discrete comparison: are we converging on the discrete sims?

%% run iterative sims to get steady-state mito distribution
options=struct();
options.nstep = 1e6;
options.Km = 0.1;
options.c0 = 0.1;
options.kw = 1;
options.L = 500;
options.ks=100;
options.kg=1;
options.dodisplay=0;
options.showevery=100;
options.nmito=75;

options.delt = 1e-6;
options.nstep = 1e7;
options.gpts = 100;
[gluc_itr,Tmito_itr,Smito_itr,Smito_int_itr,normdtg,gluc_init,opt,xpos,lmdh,ftc] = runiterativesims(options);

%%
% sample initial mitochondria positions from Tmito
Tmito = Tmito_itr;
dx = xpos(2)-xpos(1);
Nfact = trapz(Tmito_itr)*dx;
Tmito = Tmito/Nfact;

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

startpos = max(startpos,0.5/500);
startpos = min(startpos,options.L-0.5/500);

% compare histograms
[freq,bins] = hist(startpos,10)
plot(bins,freq)
db = bins(2)-bins(1);
plot(bins,freq/nsamp/db)
hold all
plot(xpos,Tmito)
hold off


%% Run discrete simulations starting with iterative solution (for comparison with analytical results)
options=struct();
options.nstep = 1e6;
options.Km = 0.1;
options.c0 = 0.1;
options.kw = 1;
options.L = 500;
options.ks=100;
options.kg=1;
options.dodisplay=0;
options.showevery=500;
options.nmito=75;
options.startgluc = gluc_itr;
options.glucfix = 1;

options.delt = 0.05;
options.gpts = 100;
nitr = 100;

clear varmito gluc_dis mitopos_dis
for j = 1:1:nitr
    
    % --- sample initial mito distribution from iterative results -----
    clear startpos
    nsamp = options.nmito;
    for mc = 1:nsamp
        u = rand();
        startpos(mc) = interp1(cTmito,xposdistrib,u);
    end
    
    startpos = max(startpos,0.5/500);
    startpos = min(startpos,options.L-0.5/500);
    options.startpos = startpos'*options.L; % dimensional starting positions
    % ----------------------------------------------------------------
    
    [gluc, mitopos, mitostate, opt] = rundiscretesims(options);
    varmito(j) = var(mitopos) ; %variance in mitochondria position distribution;
    gluc_dis(j,:) = gluc;
    mitopos_dis(j,:) = mitopos;
    
    [j varmito(j)]
    
    %save the workspace
    formatOut = 'yyyymmdd';
    date = datestr(datetime('today'),formatOut);
    %save workspace with today's date'
    filename = strcat('/home/ekoslover/proj/mitochondriaMotility/results/workspace_',date,'discretesims_100itr_1e6steps_startmitoiter_fixgluc');
    save (filename);
end