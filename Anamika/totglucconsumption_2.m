
% --------------------
%% run iterative sims to get total glucose consumption vs ks
options=struct();
options.nstep = 1e6;
options.Km = 0.1;
options.c0 = 1;
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
c0list = [0.1 0.5 1 10];
nc0 = 4;

allgluc = zeros(options.gpts,nk);
allTmito = allgluc;
consumption = zeros(1,nk,nc0);
for c0 = 1:1:nc0
    options.c0 = c0list(c0);
    for kc= 1:nk
        options.ks = kslist(kc);
        [gluc_itr,Tmito_itr,Smito_itr,Smito_int_itr,normdtg,gluc_init,opt,xpos,lmdh,ftc] = runiterativesims(options);
        allgluc(:,kc,c0) = gluc_itr;
        allTmito(:,kc,c0) = Tmito_itr;
        gluc = gluc_itr * options.Km;
        
        % get total glucose consumption
        %integrand = options.kg*options.Km*options.nmito/options.L*gluc_itr.*Tmito_itr./(gluc_itr+options.Km);
        integrand = options.kg*options.Km*options.nmito/options.L*gluc.*Tmito_itr./(gluc+options.Km);
        consumption(kc,c0) = sum((integrand(2:end)+integrand(1:end-1))/2)*(xpos(2)-xpos(1));
        [kc consumption(kc,c0)]
        %     uniform_int = options.kg*options.Km*options.nmito/options.L*gluc./(gluc+options.Km);
        %     uniconsumption(kc) = sum((uniform_int(2:end)+uniform_int(1:end-1))/2)*(xpos(2)-xpos(1));
    end
end

%%
% convert to # of gluc per second per mitochondrion
scl = 1e-3*6e23/1000*1e-12*(pi*2^2*500)/75;
for c0 = 1:1:nc0
    
    loglog(kslist,consumption(:,c0)*scl,'-')
    ylim([50000,1e6]);
    hold on;
    
end
% compare with uniform consumption