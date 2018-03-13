% Run discrete simulations starting with iterative solution (for comparison with analytical results)
options=struct();
options.nstep = 1e6;
options.Km = 0.1;
options.c0 = 0.1;
options.kw = 1;
options.L = 500;
options.ks=100;
options.kg=1;
options.dodisplay=1;
options.showevery=500;
options.nmito=75;
%options.startgluc = gluc_itr;

options.delt = 0.05;
options.gpts = 100;
nitr = 100;

% --- sample initial mito distribution from iterative results -----
%     clear startpos
%     nsamp = options.nmito;
%     for mc = 1:nsamp
%         u = rand();
%         startpos(mc) = interp1(cTmito,xposdistrib,u);
%     end
%
%     startpos = max(startpos,0.5/500);
%     startpos = min(startpos,options.L-0.5/500);
%     options.startpos = startpos'*options.L; % dimensional starting positions
% ----------------------------------------------------------------

[gluc, mitopos, mitostate, opt] = discretesims_animation(options);
varmito = var(mitopos) ; %variance in mitochondria position distribution;
gluc_dis = gluc;
mitopos_dis = mitopos;
varmetric = (6*var(mitopos)) - 0.5; %Lh = 1

