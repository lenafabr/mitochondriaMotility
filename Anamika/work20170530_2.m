
%% options for testing
%modified basicwrapper for varying c0
options = struct();
options.nmito = 14*5;

%c0 = 0.01

options.L = 500;
options.D = 140;
options.kg = 2;
%options.Km = 0.1/c0;
%options.startpos = 50;
options.pstartwalk = 1;
options.nstep = 1e4;
options.restart = 1;


options.showevery = 100;
options.dodisplay = 0; 


%options.ks = (1/4.8*1e-6)*c0*(10^-3*6e23/1000/1e12*4^2);
options.kw = 1 * 0.1;
options.delt=5e-2;

%options.startgluc = gluc;
%options.startgluc(end) = 0.1;
nc0 = 110;
% mitorec_c0 = zeros(70,nc0);
nitr = 50;
for i = 1:1:nc0
    options.c0 = i*0.01;
    c0 = i * 0.01;
%     c0 = 1;
%     options.c0 = 1;
    options.Km = 0.1/c0;
    options.ks = 10*(1/4.8*1e-6)*c0*(10^-3*6e23/1000/1e12*4^2);
    mitopos_total = [];
    for j = 1:1:nitr
        %average over the distribution, over nitr iterations
        [gluc,mitopos,mitostate,opt] = runmitosim_michaelis2(options);
        mitopos_total = [mitopos_total;mitopos];
%         mitodis = mitodis + countnos(options.L,mitopos);
    end
    mitopos_c0(:,i) = mitopos_total;
    i
%     mitodis_c0(:,i) = mitodis/nitr; %gives the average # of mito, in l belongs to [0,L]
    xpos = linspace(0,500,opt.gpts)';
end
