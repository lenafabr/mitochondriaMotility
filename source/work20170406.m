%% basic calculation: distribution with a finite first-order sink
% in symmetric geometry with fixed edge concentration


k = 1000;
D = 1;
l = 0.2;
L = 0.8;
c0 = 1;

% total rxn rate within sink

rr = 2*k*sqrt(D/k)*tanh(sqrt(k/D)*l)

% concentration at sink edge
cL = 2*D*c0/(2*k*sqrt(D/k)*tanh(sqrt(k/D)*l)*L + 2*D)

%% Numerical simulations

dx = 1e-2;
xvals = -1:dx:1;
sinkind = find(abs(xvals)<=0.2);
sinkwidth = 0.2;
k=10000;
D=1;

restart = 1;
if (restart)
    cvals = ones(size(xvals))*1;
    %cvals = abs(xvals);
    d2c = zeros(size(cvals));
    cvals0 = cvals;
end
nstep = 100000;
dt = 2e-5;
tvals = 0:dt:nstep*dt;
allcvals = zeros(size(cvals,2),nstep);
cvals0 = cvals;
for tc = 1:nstep
    d2c = zeros(size(cvals));
    d2c(2:end-1) = (cvals(3:end)+cvals(1:end-2) - 2*cvals(2:end-1))/dx^2;

    dtc = D*d2c;
    dtc(sinkind) = dtc(sinkind)-k*cvals(sinkind);
    %dtc = dtc - k *cvals.*exp(-xvals.^2/sinkwidth);
    dtc(1) = 0; dtc(end) = 0;
    cvals = cvals+dtc*dt;
    %allcvals(:,tc) = cvals';
    
    plot(xvals,cvals,'.-')%,xvals,dtc)
           ylim([0,1])
    if (mod(tc,10) == 0)
       % disp([tc,cvals(1201)])

        drawnow
    end
    %hold all
end
hold off

%%