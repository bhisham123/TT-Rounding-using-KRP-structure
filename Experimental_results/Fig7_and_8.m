close all
clear all

load('Cookies_NewData_results.mat')

%% statistics
med_sumTime_Normal       = mean(sumTimeNormalGMRES, 2);
low_sumTime_Normal       = med_sumTime_Normal - min(sumTimeNormalGMRES, [], 2);
high_sumTime_Normal       = max(sumTimeNormalGMRES, [], 2) - med_sumTime_Normal;

med_runTime_Normal       = mean(runtimeNormalGMRES, 2);
low_runTime_Normal       = med_runTime_Normal - min(runtimeNormalGMRES, [], 2);
high_runTime_Normal       = max(runtimeNormalGMRES, [], 2) - med_runTime_Normal;


med_sumTime_TTlike       = mean(sumTimeRandGMRES, 2);
low_sumTime_TTlike       = med_sumTime_TTlike - min(sumTimeRandGMRES, [], 2);
high_sumTime_TTlike       = max(sumTimeRandGMRES, [], 2) - med_sumTime_TTlike;

med_runTime_TTlike       = mean(runtimeRandGMRES, 2);
low_runTime_TTlike       = med_runTime_TTlike - min(runtimeRandGMRES, [], 2);
high_runTime_TTlike       = max(runtimeRandGMRES, [], 2) - med_runTime_TTlike;


med_sumTime_KRP       = mean(sumTimeRandKRPGMRES, 2);
low_sumTime_KRP       = med_sumTime_KRP - min(sumTimeRandKRPGMRES, [], 2);
high_sumTime_KRP       = max(sumTimeRandKRPGMRES, [], 2) - med_sumTime_KRP;

med_runTime_KRP       = mean(runtimeRandKRPGMRES, 2);
low_runTime_KRP       = med_runTime_KRP - min(runtimeRandKRPGMRES, [], 2);
high_runTime_KRP       = max(runtimeRandKRPGMRES, [], 2) - med_runTime_KRP;


%% Time and rank plot
f=figure('Position', [100, 100, 1100, 400]);

%Time plot
sp1 = subplot(1,2,1);
hold on;

sumtime_p1 = errorbar(N(1:4), med_sumTime_Normal, low_sumTime_Normal, high_sumTime_Normal, 'ko-', 'LineWidth', 2.2,'MarkerSize', 10, 'DisplayName', 'Sum-Round (Naive)');
% runtime_p1 = errorbar(N(1:4), med_runTime_Normal, low_runTime_Normal, high_runTime_Normal, '-s','Color' ,'#84563C','LineWidth', 2.2,'MarkerSize', 10, 'DisplayName', 'Total time (Naive)');
runtime_p1 = errorbar(N(1:4), med_runTime_Normal, low_runTime_Normal, high_runTime_Normal, '-.s','Color' ,'#808080','LineWidth', 2.2,'MarkerSize', 10, 'DisplayName', 'Total time (Naive)');


runtime_p2 = errorbar(N, med_runTime_KRP, low_runTime_KRP, high_runTime_KRP, '-.s','Color','#FF8080','LineWidth', 2.2,'MarkerSize', 10, 'DisplayName','Total time (Rand-Orth-KRP)');
sumtime_p2 = errorbar(N, med_sumTime_KRP, low_sumTime_KRP, high_sumTime_KRP, 'ro-','LineWidth', 2.2,'MarkerSize', 10, 'DisplayName','Sum-Round (Rand-Orth-KRP)');
% runtime_p2 = errorbar(N, med_runTime_KRP, low_runTime_KRP, high_runTime_KRP, '-s','Color','#FF8C00','LineWidth', 2.2,'MarkerSize', 10, 'DisplayName','Total time (Rand-Orth-KRP)');

runtime_p3 = errorbar(N, med_runTime_TTlike, low_runTime_TTlike, high_runTime_TTlike, '-.s','Color','#8080FF','LineWidth', 2.2,'MarkerSize', 10, 'DisplayName', 'Total time (Rand-Orth)');
sumtime_p3 = errorbar(N, med_sumTime_TTlike, low_sumTime_TTlike, high_sumTime_TTlike, 'bo-','LineWidth', 2.2,'MarkerSize', 10, 'DisplayName', 'Sum-Round (Rand-Orth)');
% runtime_p3 = errorbar(N, med_runTime_TTlike, low_runTime_TTlike, high_runTime_TTlike, '-s','Color','#7E2F8E','LineWidth', 2.2,'MarkerSize', 10, 'DisplayName', 'Total time (Rand Orth)');


set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('Number of parameter samples n_2 = \ldots = n_d (log scale)', 'FontSize', 14)
ylabel('Time (Seconds)', 'FontSize', 14)
xticks(N);
ylim([5*10^1, 10^5])
set(gca,'FontSize',14)
lgd1 = legend(sp1, [runtime_p2,sumtime_p2,runtime_p3,sumtime_p3,runtime_p1,sumtime_p1]);
set(lgd1,'Units', 'normalized','Position',[0.278181818181818 0.1225 0.185 0.24125], 'FontSize', 12);
box on
grid on;

%rank plot
subplot(1,2,2);
hold on
ranks_krp = max(ranksrandkrp{1}(1:end-1,:),[],2);
p3 = plot(ranks_krp, 'r*','markersize',10,'linewidth',2);
ranks_rand = max(ranksrand{1}(1:end-1,:),[],2);
p1 = plot(ranks_rand, 'bo','markersize',10,'linewidth',2);
ranks_det = max(ranksnormal{1}(1:end-1,:),[],2);
p2 = plot(ranks_det, 'k+','markersize',10,'linewidth',2);
for j=2:length(N)
    ranks = max(ranksrandkrp{j}(1:end-1,:),[],2);
    plot(ranks, 'r*', 'markersize',10,'linewidth',2)
end
for j=2:length(N)
    ranks = max(ranksrand{j}(1:end-1,:),[],2);
    plot(ranks, 'bo', 'markersize',10,'linewidth',2)
end
for j=2:length(N)
    ranks = max(ranksnormal{j}(1:end-1,:),[],2);
    plot(ranks, 'k+', 'markersize',10,'linewidth',2)
end
hold off

xlabel('TT-GMRES iteration number', 'FontSize', 14)
ylabel('TT-rank', 'FontSize', 14)
set(gca,'FontSize',14)
legend2 = legend([p3,p1,p2], 'TT-GMRES-Rand-Orth-KRP' , 'TT-GMRES-Rand-Orth','Naive TT-GMRES');

set(legend2,'Position',[0.729090909090909 0.12125 0.175 0.12125],'FontSize',12);
box on
grid on;
exportgraphics(f, 'Fig8.png','Resolution', 400)


%% Relative Error plot 
avg_err_normal = mean(res_err_normal,2);
avg_err_rand =  mean(res_err_rand,2);
avg_err_randKRP = mean(res_err_randKRP,2);

f = figure('Position', [100, 100, 520, 400]);
hold on 
plot(N,avg_err_randKRP(1:length(N)),'r*-', 'LineWidth',2.2,'markersize',10,'DisplayName','TT-GMRES-Rand-Orth-KRP'); 
plot(N,avg_err_rand(1:length(N)),'bo-','LineWidth',2.2,'markersize',10,'DisplayName','TT-GMRES-Rand-Orth')
plot(N,avg_err_normal(1:length(N)),'k+-','LineWidth',2.2,'markersize',10,'DisplayName','Naive TT-GMRES')

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

xlabel('Number of parameter samples n_2 = \ldots = n_d (log scale)', 'FontSize', 14)
ylabel('Relative residual error', 'FontSize', 14)
xticks(N);
set(gca,'FontSize',15)

lgd = legend('FontSize', 11.8);
set(lgd, 'Units', 'normalized');                     % Use normalized units
set(lgd, 'Position', [0.4, 0.45, 0.2, 0.1]);          % Manually center it

box on
grid on
exportgraphics(f, 'Fig7.png','Resolution', 400)



