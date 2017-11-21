% compare analytical solutions for fixed conc vs fixed permeability

xlist = linspace(-0.5,0.5);
lam=0.1;
D=140;
p=100/D;

% fixed concentration
Gconc = cosh((xlist)/lam)/cosh(1/2/lam);
%plot(xlist,Gconc)

% permeability
k=1/lam;
B =  (p*(k*cosh(k/2) + p*sinh(k/2))) / ((k^2 + p^2)*sinh(k) + 2*k*p*cosh(k));
A = (B*(k-p)*exp(-k) + p*exp(-k/2))/(k+p);
Gperm = A*exp(k*xlist) + B*exp(-k*xlist);
plot(xlist,Gconc)
hold all
plot(xlist,Gperm)
hold off

