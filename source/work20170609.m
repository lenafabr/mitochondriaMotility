% variance metric as a function of two dimless parameters: 
% lh = sqrt(D/(kg*N*del*L))
% Al = ks*c0/kw * lh
% assumes c0 is small enough that glucose consumption rates are linear
% throughout (ie: c0<<Km)
% assumes glucose distribution is similar to that expected with uniform
% mitochondrial density

varfunc  = @(lh,Al) 2*Al*(-6*lh + (1+12*lh^2)*tanh(1/2/lh)) / (1+2*Al*tanh(1/2/lh));

lhlist = logspace(-3,1,51);
Alist = logspace(-2,2,50);

varmetric = zeros(length(lhlist),length(Alist));
for lc = 1:length(lhlist)
    for ac = 1:length(Alist)
        varmetric(lc,ac) = varfunc(lhlist(lc),Alist(ac));
    end
end


pcolor(log10(Alist),log10(lhlist),varmetric)
shading flat
%xlabel('A*lh')
%ylabel('lh = sqrt(D/(kg*N*del*L))')
title('variance metric')

%%
% get relevant numbers
D=140; L=500; N = 70; del=1;
kg = 0.2;
lh = sqrt(D/(kg*N*del*L))

%%
hold all
plot(log10(Alist),log10(lh)*ones(size(Alist)),'k','LineWidth',2)
hold off

%% for a given value of A, plot versus lh
Alist = [1,10,100];
lhlist = logspace(-3,1,51);
varmetric = zeros(length(lhlist),length(Alist))
for ac = 1:length(Alist)
    A = Alist(ac);
    for lc = 1:length(lhlist)
        lh = lhlist(lc);
        varmetric(lc,ac) = varfunc(lh,A*lh);
    end
    semilogx(lhlist,varmetric(:,ac),'LineWidth',2)
    hold all
end
hold off

%xlabel('lh')
%ylabel('variance metric')
legend('A=1','A=10','A=100')
set(gca,'FontSize',14)

%% for a given value of lh, plot vs A
Alist = logspace(-2,2);;
lhlist = [0.01,0.05,0.1,0.2,1];
varmetric = zeros(length(lhlist),length(Alist))
for lc = 1:length(lhlist)
    for ac = 1:length(Alist)
        varmetric(lc,ac) = varfunc(lhlist(lc),Alist(ac)*lhlist(lc));
    end
    semilogx(Alist,varmetric(lc,:))
    hold all
    strlabel{lc} = sprintf('lh=%0.2f',lhlist(lc));
end
hold off

xlabel('A')
ylabel('variance metric')
legend(strlabel)

%% sim data from Anamika
load('results/workspace_6_12_50itr_10fold.mat')
