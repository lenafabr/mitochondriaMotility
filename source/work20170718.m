load('workspace_07_10_1')
var_metric = var_metric*6/opt.L^2 - 1/2;
%% surface plot results
pcolor(log10(A_var*opt.kw/opt.ks),log10(Lambda_hat),var_metric); shading flat

%% horizontal slices (versus A of versus c0)
plot(log10(A),var_metric(10:10:100,:)')
xlabel('A')
ylabel('variance')
%%
plot(A,var_metric(10:10:100,:)')

%% vertical slices (versus lh)
plot(log10(lmdh),var_metric(:,10:10:100))
xlabel('lambda hat')
ylabel('variance')

%% get optimal concentrations
for lc = 1:length(lmdh)
    [varmax(lc),ind] = max(var_metric(lc,:));
    Aopt(lc) = A(ind)*opt.kw/opt.ks;
end

loglog(lmdh,Aopt,lmdh,4e-3./lmdh.^2)
xlabel('lambda hat')
ylabel('optimal c0/km')