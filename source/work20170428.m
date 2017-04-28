for c = 1:length(mitopos)
    if (mitopos(c)>100 & mitopos(c)<200)
        temppos(c) = rand()*100;
    elseif (mitopos(c)>200 & mitopos(c)<300)
        temppos(c) = rand()*100+300;
    else
        temppos(c) = mitopos(c);
    end
end

%%
glucm = interp1(xpos,gluc,temppos);

corr(temppos',glucm')

%%
[nmito,xbin] = hist(temppos,40);
plot(xbin,nmito)

glucbin = interp1(xpos,gluc,xbin);

corr(glucbin',nmito')