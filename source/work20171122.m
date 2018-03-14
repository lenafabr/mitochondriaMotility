load('workspace_20171124p_kg_101_102.mat')
%% view results for constant permeability, vs kg and P

enrichind = find(xpos<51);

for kc = 1:nkg
    for pc = 1:nP
        Tmito = Tmito_all(:,kc,pc);
        mitomean(kc,pc) = trapz(Tmito.*xpos)*(xpos(2)-xpos(1));
        mitomean(kc,pc) = (mitomean(kc,pc)-250)/250;
               
        normfact = 1/500*xpos(enrichind(end));
        enrichment(kc,pc) = trapz(Tmito(enrichind))*(xpos(2)-xpos(1))/normfact;
    end
end

%%
pcolor(log10(Plist),log10(kglist(1:nkg)),mitomean)
xlabel('permeability')
ylabel('kg')
shading flat
colormap jet

%%
pcolor(log10(Plist),log10(kglist(1:nkg)),enrichment)
xlabel('permeability')
ylabel('kg')
shading flat
hold all
plot(log10(Plist),log10(Plist*500/70*25),'k--')
hold off
colormap jet

% plot using our estimates of kg and P
hold all
plot(log10(0.05),log10(1.05),'k.','MarkerSize',20)
hold off
%% look at tmito distribution
plot(xpos,Tmito_all(:,1,66),xpos,1/500*ones(size(xpos)),'--')

%% compare distributions across boundary
cmap = colormapinterp([0,0,1;1,0,0],43);
for c = 1:43
plot(xpos,Tmito_all(:,39,c),'Color',cmap(c,:))
hold all
end

hold off