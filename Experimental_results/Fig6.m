close all
clear all

load('Cookies_OldData_results.mat')

f = figure(1);
f.Position(1:2) = [0, 0];
f.Position(3:4) = [1350, 350];
X = categorical(N);


tiledlayout(1,3);

nexttile;
bar(X, [mean(sumTimeRandKRPGMRES,2), mean(otherTimeRandKRPGMRES,2)], 'stacked')
title('Randomized TT-GMRES-KRP')
xlabel('Number of parameter samples', 'FontSize', 18)
ylabel('Average Runtime (s)', 'FontSize', 18)
legend('RandRoundSum', 'Other', 'location', 'northwest')
set(gca,'FontSize',16)

nexttile;
bar(X, [mean(sumTimeRandGMRES,2), mean(otherTimeRandGMRES,2)], 'stacked')
title('Randomized TT-GMRES')
xlabel('Number of parameter samples', 'FontSize', 18)
ylabel('Average Runtime (s)', 'FontSize', 18)
legend('RandRoundSum', 'Other', 'location', 'northwest')
set(gca,'FontSize',16)

nexttile;
bar(X, [mean(sumTimeNormalGMRES,2), mean(otherTimeNormalGMRES,2)], 'stacked')
title('Naive TT-GMRES')
xlabel('Number of parameter samples', 'FontSize', 18)
ylabel('Average Runtime (s)', 'FontSize', 18)
legend('TT-Sum + Round', 'Other', 'location', 'northwest')
set(gca,'FontSize',16)



f=figure('Position', [100, 100, 1100, 400]);
subplot(1,2,1);


%Randomized with KRP Structure

speedup_krp = median(sumTimeNormalGMRES, 2) ./ median(sumTimeRandKRPGMRES, 2);
neg_krp = speedup_krp - min(sumTimeNormalGMRES,[],2)./ max(sumTimeRandKRPGMRES,[],2);
pos_krp = max(sumTimeNormalGMRES,[],2)./ min(sumTimeRandKRPGMRES,[],2) - speedup_krp;
errorbar(N, speedup_krp, neg_krp, pos_krp,'r-o','markersize',10,'linewidth',2);
hold on
speedup_krp = median(runtimeNormalGMRES, 2) ./ median(runtimeRandKRPGMRES, 2);
neg_krp = speedup_krp - min(runtimeNormalGMRES,[],2)./ max(runtimeRandKRPGMRES,[],2);
pos_krp = max(runtimeNormalGMRES,[],2)./ min(runtimeRandKRPGMRES,[],2) - speedup_krp;
% errorbar(N, speedup_krp, neg_krp, pos_krp,'-s','Color','#FF8C00','markersize',10,'linewidth',2);
errorbar(N, speedup_krp, neg_krp, pos_krp,'-.s','Color','#FF8080','markersize',10,'linewidth',2);





%Randomized with TT Structure
speedup = median(sumTimeNormalGMRES, 2) ./ median(sumTimeRandGMRES, 2);
neg = speedup - min(sumTimeNormalGMRES,[],2)./ max(sumTimeRandGMRES,[],2);
pos = max(sumTimeNormalGMRES,[],2)./ min(sumTimeRandGMRES,[],2) - speedup;
errorbar(N, speedup, neg, pos,'b-o','markersize',10,'linewidth',2);
hold on 
speedup = median(runtimeNormalGMRES, 2) ./ median(runtimeRandGMRES, 2);
neg = speedup - min(runtimeNormalGMRES,[],2)./ max(runtimeRandGMRES,[],2);
pos = max(runtimeNormalGMRES,[],2)./ min(runtimeRandGMRES,[],2) - speedup;
% errorbar(N, speedup, neg, pos,'-s','Color','#7E2F8E','markersize',10,'linewidth',2);  %'Color','#008000' 
errorbar(N, speedup, neg, pos,'-.s','Color','#8080FF','markersize',10,'linewidth',2)


set(gca, 'XScale', 'log')
xlabel('Number of parameter samples n_2 = \ldots = n_d (log scale)', 'FontSize', 14)
ylabel('Speedup', 'FontSize', 14)
xticks(N);
set(gca,'FontSize',14)
lgd1=legend('Sum-Round (Rand-Orth-KRP)','Total speedup (Rand-Orth-KRP)', 'Sum-Round (Rand-Orth)','Total speedup (Rand-Orth)');
set(lgd1,'Position',[0.270548295454545 0.1208125 0.193636363636364 0.15875],'FontSize',12);
box on
grid on


subplot(1,2,2);
ranks_krp = max(ranksrandkrp{1}(1:end-1,:),[],2);
p3 = plot(ranks_krp, 'r*','markersize',10,'linewidth',2); hold on
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
lgd2 = legend([p3,p1,p2], 'TT-GMRES-Rand-Orth-KRP' , 'TT-GMRES-Rand-Orth','Naive TT-GMRES');
set(lgd2,'Position',[0.73 0.12 0.175 0.12125],'FontSize',12);
box on
grid on

exportgraphics(f, 'Fig6.png','Resolution', 400)


avg_err_normal = mean(res_err_normal,2);
avg_err_rand =  mean(res_err_rand,2);
avg_err_randKRP = mean(res_err_randKRP,2);

f = figure('Position', [100, 100, 520, 400]);
hold on 
plot(N,avg_err_randKRP(1:length(N)),'r*-', 'LineWidth',2,'markersize',10,'DisplayName','TT-GMRES-Rand-Orth-KRP'); 
plot(N,avg_err_rand(1:length(N)),'bo-','LineWidth',2,'markersize',10,'DisplayName','TT-GMRES-Rand-Orth')
plot(N,avg_err_normal(1:length(N)),'k+-','LineWidth',2,'markersize',10,'DisplayName','Naive TT-GMRES')

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

xlabel('Number of parameter samples n_2 = \ldots = n_d (log scale)', 'FontSize', 14)
ylabel('Relative Error', 'FontSize', 14)
xticks(N);
set(gca,'FontSize',14)

legend('location', 'best','FontSize', 12.5)

box on
grid on


