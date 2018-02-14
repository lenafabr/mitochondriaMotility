%% run iterative sims to get total glucose consumption vs ks
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
%%
nk=30;
kslist = logspace(0,4,nk);

allgluc = zeros(options.gpts,nk);
allTmito = allgluc;
consumption = zeros(1,nk);
for kc= 1:nk
    kc
    options.ks = kslist(kc);
    [gluc_itr,Tmito_itr,Smito_itr,Smito_int_itr,normdtg,gluc_init,opt,xpos,lmdh,ftc] = runiterativesims(options);
    allgluc(:,kc) = gluc_itr;
    allTmito(:,kc) = Tmito_itr;
    
    % get total glucose consumption
    integrand = options.kg*options.Km*options.nmito/options.L*gluc_itr.*Tmito_itr./(gluc_itr+options.Km);
    consumption(kc) = sum((integrand(2:end)+integrand(1:end-1))/2)*(xpos(2)-xpos(1));
    [kc consumption(kc)]
end

%%
% convert to # of gluc per second per mitochondrion
scl = opt.Km*1e-3*6e23/1000*1e-12*(pi*2^2*500)/75;
semilogx(kslist*opt.Km,consumption*scl,'-')