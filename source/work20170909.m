%load('./Anamika/workspace_20170901uniform.mat')
%load('./Anamika/workspace_20170901edges.mat')
load('../Anamika/workspace_20170902continuous_sol.mat')
%%
% average variance metric
mean(varmito)*6/options.L^2-0.5
%std(varmito)*6/options.L^2 / sqrt(100)

%% overall variance metric
allpos = reshape(mitopos_dis,1e4,1);
varall = var(allpos)
varall*6/options.L^2 - 0.5