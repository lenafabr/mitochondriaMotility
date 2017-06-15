% Michaelis-menten glucose consumption kinetics
%% options for testing
options = struct();
options.nmito = 14*5;

c0 = 0.4

options.L = 500;
options.D = 140;
options.kg = 0.2;
options.Km = 0.1/c0;
%options.startpos = 50;
options.pstartwalk = 1;
options.nstep = 1e4;

options.showevery = 100;

options.cend=1;
options.ks = (1/4.8*1e-6)*c0*(10^-3*6e23/1000/1e12*4^2);
options.kw = 0.1;
options.delt=5e-2;
%
%options.startgluc = gluc;
%options.startgluc(end) = 0.1;
%options.fixgluc = gluc;

[gluc,mitopos,mitostate,opt] = runmitosim_michaelis2(options)
xpos = linspace(0,500,opt.gpts)';

%% upper limit based on fixed glucose profile, linear kinetics

lh = sqrt(options.D/(options.kg*options.L*options.nmito))
A = options.ks/options.kw
varfunc  = @(lh,Al) 2*Al*(-6*lh + (1+12*lh^2)*tanh(1/2/lh)) / (1+2*Al*tanh(1/2/lh));
varfunc(lh,A*lh)
%%
save('~/grants/NSFcareer/workfig/mitochondriasim.mat')
%% variance metric
(var(mitopos)-options.L^2/12)/(options.L^2/6)
%% equilibrium constant for binding
kstop = options.ks*options.Km*gluc./(options.Km+gluc);
kstart = options.kw;
Keq = kstop./kstart;

%% plot simulation results
plot(xpos,gluc,'LineWidth',2)
hold all
walkind = find(mitostate);
stopind = find(~mitostate);
plot(mitopos(walkind),0.5*ones(size(walkind)),'o','Color',[0,0.5,0],'MarkerSize',5,'LineWidth',2)
plot(mitopos(stopind),0.5*ones(size(stopind)),'o','Color','r','MarkerSize',5,'LineWidth',2)
hold off
axis square
%set(gca,'Position',[0.05 0.05 0.9 0.9],'XTickLabel',[],'YTickLabel',[],'LineWidth',1)

%% plot analytical mitochondria steady-state distribution, given fixed G(x)
kstop = options.ks*options.Km*gluc./(options.Km+gluc);
kstart = options.kw;
Keq = kstop./kstart;
Keqint = sum((Keq(2:end)+Keq(1:end-1))/2)*(xpos(2)-xpos(1));
L = options.L
mdistrib = 1/(L+Keqint)*(Keq+1);
plot(xpos,mdistrib*L)

%% plot analytical glucose solution for uniform mitochondria, linear kinetics
kg = options.kg; L = options.L; D = options.D;
f = options.nmito/L;
func = @(x) cosh(sqrt(kg*f/D)*(x-L/2))/cosh(sqrt(kg*f/D)*L/2);
hold all
plot(xpos,func(xpos))
hold off
%% histogram
hist(mitopos,9)
axis square
%set(gca,'Position',[0.05 0.05 0.9 0.9],'XTickLabel',[],'YTickLabel',[])
