%% attempt to reproduce old results using new (paper) parameters
% discrete sims with gpts = 500
%only 1 iteration, just to check

options.nstep = 1e8; 
options.Km = 0.1;
options.c0 = 0.1;
options.kw = 1;
options.L = 500;
options.ks=100;
options.kg=1;
options.dodisplay=0;
options.showevery=100;
options.nmito=75;
options.dttol = 1e-3;

options.delt = 1e-3;
options.gpts = 500;
%nitr = 100;
nitr = 1;
clear varmito gluc_dis mitopos_dis
for j = 1:1:nitr
    [gluc, mitopos, mitostate, opt] = rundiscretesims(options);
    varmito(j) = var(mitopos) ; %variance in mitochondria position distribution;
    gluc_dis(j,:) = gluc;
    mitopos_dis(j,:) = mitopos;
    
    [j varmito(j)]
end
%save the workspace
formatOut = 'yyyymmdd';
date = datestr(datetime('today'),formatOut);
%save workspace with today's date'
filename = strcat('workspace_',date,'discretesims_100itr_gpts500');
save (filename);

%% concatenate mitopos_dis to one long array
mitodis = reshape(mitopos_dis,[1,size(mitopos_dis,1)*size(mitopos_dis,2)]);
hist_mito = histogram(mitodis,20,'Normalization','pdf');
%hold on
%histogram(mitodis,100,'Normalization','pdf');
%hold off;




%% then generate plot from iterative sims
options.Km = 0.1;
options.c0 = 0.1;
options.kw = 1;
options.L = 500;
options.ks=100;
options.kg=1/10;
options.dodisplay=0;
options.showevery=100;
options.nmito=750;

options.delt = 1e-6;
options.nstep = 1e7;
options.gpts = 500;
[gluc_itr,Tmito_itr,Smito_itr,Smito_int_itr,normdtg,gluc_init,opt,xpos,lmdh,ftc] = runiterativesims(options);

dx = xpos(2)-xpos(1);
Nfact = trapz(Tmito_itr)*dx;
hold all
plot(xpos,Tmito_itr/Nfact)
hold off