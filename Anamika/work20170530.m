%script to run iterations for the same parameters 
%to get averaged-out mitopos 
nitr = 10;
mitorec = zeros(70,1);
mitovar = 0;
for i = 1:1:nitr
    basicrunwrapper
    mitorec = mitorec + mitopos;
    mitovar = mitovar + var(mitopos);
end
mitorec = mitorec /nitr;

