clc
clear all 
close all

load('Matern_results.mat')

mean_errLR = mean(errorsLR,2);
mean_errRandLR = mean(errorsRandLR,2);
mean_errRandKtr = mean(errorsRandKRP,2);

mean_timeLR = mean(timeLR(:,2:end),2);
mean_timeRandLR = mean(timeRandLR(:,2:end),2);
mean_timeRandKtr = mean(timeRandKRP(:,2:end),2);

mean_timeRandLR_fix = mean(timeRandLR_fix(:,2:end),2);
mean_timeRandOrth_fix = mean(timeRandOrth_fix(:,2:end),2);
mean_timeRandKtr_fix = mean(timeRandKRP_fix(:,2:end),2);

mean_speedupLR = mean_timeLR./mean_timeLR;
mean_speedupRandLR = mean_timeLR./mean_timeRandLR;
mean_speedupRandKtr = mean_timeLR./mean_timeRandKtr;

mean_speedupRandLR_fix = mean_timeLR./mean_timeRandLR_fix;
mean_speedupRandOrth_fix = mean_timeLR./mean_timeRandOrth_fix;
mean_speedupRandKtr_fix = mean_timeLR./mean_timeRandKtr_fix;


mean_compLR = mean(compLR,2);
mean_compRandLR = mean(compRandLR,2);
mean_compRandKRP = mean(compRandKRP,2);


f = figure('Position', [0, 0, 1000, 4000]);
% Compute statistics
med_LR       = mean(errorsLR, 2);
min_LR       = min(errorsLR, [], 2);
max_LR       = max(errorsLR, [], 2);

med_RandLR   = mean(errorsRandLR, 2);
min_RandLR   = min(errorsRandLR, [], 2);
max_RandLR   = max(errorsRandLR, [], 2);

med_RandKtr  = mean(errorsRandKRP, 2);
min_RandKtr  = min(errorsRandKRP, [], 2);
max_RandKtr  = max(errorsRandKRP, [], 2);

% Compute error bar lengths (top and bottom)
err_low_LR     = med_LR - min_LR;
err_high_LR    = max_LR - med_LR;

err_low_RandLR  = med_RandLR - min_RandLR;
err_high_RandLR = max_RandLR - med_RandLR;

err_low_RandKtr  = med_RandKtr - min_RandKtr;
err_high_RandKtr = max_RandKtr - med_RandKtr;

% X-axis positions
x = 1:size(errorsLR, 1);  % Assuming same number of rows for all


sp1 = subplot(3,2,1,'Parent',f);
hold on;

p1 = errorbar(5*x, med_LR, err_low_LR, err_high_LR, 'ko', 'LineWidth', 2,'MarkerSize', 10, 'DisplayName', 'TT-Rounding');
p2 = errorbar(5*x, med_RandKtr, err_low_RandKtr, err_high_RandKtr, 'r*','LineWidth', 2,'MarkerSize', 10, 'DisplayName','TT-KRP (Adap)');
p3 = errorbar(5*x, med_RandLR, err_low_RandLR, err_high_RandLR, 'b+','LineWidth', 2,'MarkerSize', 10, 'DisplayName', 'TT-Orth-Rand (Adap)');

xticks(5*x)
xticklabels({'1e-1','1e-2', '1e-3', '1e-4', '1e-5','1e-6','1e-7'})
xlim([3,37])
ylim([1e-12, 2e-1])
yticks([1e-11,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1])
xlabel('Tolerance', 'FontSize', 14)
ylabel('Relative Error', 'FontSize', 14)
set(sp1, 'FontSize', 14)
set(sp1, 'YScale', 'log');
lgd1 = legend(sp1, [p2, p3, p1]);
box on
grid on;


sp2 = subplot(3,2,3,'Parent',f);
hold on
plot(5*(1:length(mean_speedupRandKtr)), mean_speedupRandKtr,   'r*-','markersize',10.5,'linewidth',2,'DisplayName', 'TT-KRP (Adap)'); 
plot(5*(1:length(mean_speedupRandKtr_fix)), mean_speedupRandKtr_fix,   's--','Color','#FF8C00','markersize',10.5,'linewidth',2,'DisplayName', 'TT-KRP (Fix)'); 
plot(5*(1:length(mean_speedupRandOrth_fix)),  mean_speedupRandOrth_fix,'x--','Color','#7E2F8E','markersize',10.5,'linewidth',2,'DisplayName', 'TT-Rand-Orth')
plot(5*(1:length(mean_speedupRandLR)),  mean_speedupRandLR,'b+-','markersize',10.5,'linewidth',2,'DisplayName', 'TT-Orth-Rand (Adap)')
plot(5*(1:length(mean_speedupRandLR_fix)),  mean_speedupRandLR_fix,'bd--','markersize',10.5,'linewidth',2,'DisplayName', 'TT-Orth-Rand (Fix)')
plot(5*(1:length(mean_speedupLR)), mean_speedupLR,'ko-','markersize',10.5,'linewidth',2,'DisplayName', 'TT-Rounding')

xticks(5*(1:7))
xticklabels({'1e-1','1e-2', '1e-3', '1e-4', '1e-5','1e-6','1e-7'})
ylim([0 20])
ylabel('Speedup', 'FontSize', 14)
xlabel('Tolerance', 'FontSize', 14)
set(sp2, 'FontSize', 14)
xlim([3,37])
lgd2 = legend(sp2);
box on
grid on;
hold off

sp3 = subplot(3,2,5,'Parent',f);
plot(5*(1:length(mean_compRandKRP)), mean_compRandKRP,   'r*-','markersize',10.5,'linewidth',2,'DisplayName', 'TT-KRP (Adap)'); hold on
plot(5*(1:length(mean_compRandLR)),  mean_compRandLR,'b+-','markersize',10.5,'linewidth',2,'DisplayName', 'TT-Orth-Rand (Adap)')
plot(5*(1:length(mean_compLR)), mean_compLR,'ko-','markersize',10.5,'linewidth',2,'DisplayName', 'TT-Rounding')

xticks(5*(1:7))
xticklabels({'1e-1','1e-2', '1e-3', '1e-4', '1e-5','1e-6','1e-7'})
yticks([1e0, 1e1,1e2,1e3,1e4,1e5])
ylabel('Compression', 'FontSize', 14)
xlabel('Tolerance', 'FontSize', 14)
set(sp3, 'FontSize', 14)
set(sp3, 'YScale', 'log')
xlim([3,37])
lgd3 = legend(sp3);
box on;
grid on;
hold off;



%%% plot for Randomized rounding + Compression pass at end
mean_timeRandLR_round = mean(timeRandLR(:,2:end),2);
mean_timeRandKtr_round = mean(timeRandKRP(:,2:end)+timeRandKRP_round(:,2:end),2);

mean_speedupRandLR_round = mean_timeLR./(mean_timeRandLR_round);
mean_speedupRandKtr_round = mean_timeLR./(mean_timeRandKtr_round);

mean_compRandLR_round = mean(compRandLR_round,2);
mean_compRandKRP_round = mean(compRandKRP_round,2);

% Compute statistics
med_LR       = mean(errorsLR, 2);
min_LR       = min(errorsLR, [], 2);
max_LR       = max(errorsLR, [], 2);

med_RandLR   = mean(errorsRandLR, 2);
min_RandLR   = min(errorsRandLR, [], 2);
max_RandLR   = max(errorsRandLR, [], 2);

med_RandLR_round   = mean(errorsRandLR_round, 2);
min_RandLR_round   = min(errorsRandLR_round, [], 2);
max_RandLR_round   = max(errorsRandLR_round, [], 2);


med_RandKtr  = mean(errorsRandKRP_round, 2);
min_RandKtr  = min(errorsRandKRP_round, [], 2);
max_RandKtr  = max(errorsRandKRP_round, [], 2);

% Compute error bar lengths (top and bottom)
err_low_LR     = med_LR - min_LR;
err_high_LR    = max_LR - med_LR;

err_low_RandLR  = med_RandLR - min_RandLR;
err_high_RandLR = max_RandLR - med_RandLR;

err_low_RandLR_round  = med_RandLR_round - min_RandLR_round;
err_high_RandLR_round = max_RandLR_round - med_RandLR_round;

err_low_RandKtr  = med_RandKtr - min_RandKtr;
err_high_RandKtr = max_RandKtr - med_RandKtr;

% X-axis positions
x = 1:size(errorsLR, 1);  % Assuming same number of rows for all

sp4 = subplot(3,2,2,'Parent',f);
hold on;
p1 = errorbar(5*x, med_LR, err_low_LR, err_high_LR, 'ko', 'LineWidth', 2,'MarkerSize', 10, 'DisplayName', 'TT-Rounding');
p2= errorbar(5*x, med_RandKtr, err_low_RandKtr, err_high_RandKtr, 's','Color','#DC143C','LineWidth', 2,'MarkerSize', 10,   'DisplayName','TT-KRP (Adap-R)');
p3 = errorbar(5*x, med_RandLR, err_low_RandLR, err_high_RandLR, 'b+','LineWidth', 2,'MarkerSize', 10,  'DisplayName', 'TT-Orth-Rand (Adap)');

xticks(5*x)
xticklabels({'1e-1','1e-2', '1e-3', '1e-4', '1e-5','1e-6','1e-7'})
xlim([3,37])
ylim([1e-12, 2e-1])
yticks([1e-11,1e-10, 1e-9, 1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1])
xlabel('Tolerance', 'FontSize', 14)
ylabel('Relative Error', 'FontSize', 14)
set(sp4, 'FontSize', 14)
set(sp4, 'YScale', 'log');
lgd4 = legend(sp4,[p2,p3,p1]);
grid on;
box on;
hold off

sp5 = subplot(3,2,4,'Parent',f);
hold on
plot(5*(1:length(mean_speedupRandKtr)), mean_speedupRandKtr,   'r*-','markersize',10.5,'linewidth',2,'DisplayName', 'TT-KRP (Adap)'); hold on
plot(5*(1:length(mean_speedupRandKtr_round)), mean_speedupRandKtr_round,   's--','Color','#DC143C','markersize',10.5,'linewidth',2,'DisplayName', 'TT-KRP (Adap-R)'); hold on
plot(5*(1:length(mean_speedupRandLR)),  mean_speedupRandLR,'b+-','markersize',10.5,'linewidth',2,'DisplayName', 'TT-Orth-Rand (Adap)')
plot(5*(1:length(mean_speedupLR)), mean_speedupLR,'ko-','markersize',10.5,'linewidth',2,'DisplayName', 'TT-Rounding')

xticks(5*(1:7))
xticklabels({'1e-1','1e-2', '1e-3', '1e-4', '1e-5','1e-6','1e-7'})
ylabel('Speedup', 'FontSize', 14)
xlabel('Tolerance', 'FontSize', 14)
set(sp5, 'FontSize', 14)
xlim([3,37])
ylim([0 20])
lgd5=legend(sp5);
grid on
box on;
hold off

sp6 = subplot(3,2,6,'Parent',f);
p1 = plot(5*(1:length(mean_compLR)), mean_compLR,'ko-','markersize',10.5,'linewidth',2,'DisplayName', 'TT-Rounding');hold on
p2 = plot(5*(1:length(mean_compRandKRP_round)), mean_compRandKRP_round,   's--','Color','#DC143C','markersize',10.5,'linewidth',2,'DisplayName', 'TT-KRP (Adap-R)'); 
p3 = plot(5*(1:length(mean_compRandLR)),  mean_compRandLR,'b+-','markersize',10.5,'linewidth',2,'DisplayName', 'TT-Orth-Rand (Adap)');

xticks(5*(1:7))
xticklabels({'1e-1','1e-2', '1e-3', '1e-4', '1e-5','1e-6','1e-7'})
yticks([1e0, 1e1,1e2,1e3,1e4,1e5])
ylim([1e0,1e5])

ylabel('Compression', 'FontSize', 14)
xlabel('Tolerance', 'FontSize', 14)
set(sp6, 'FontSize', 14)
set(sp6, 'YScale', 'log')
xlim([3,37])
lgd6 = legend(sp6,[p2,p3,p1]);
grid on
box on 
hold off

drawnow;
set(lgd1,'Units', 'normalized','Position',[0.1335 0.710387399463807 0.159 0.0669999999999999], 'FontSize', 12);
set(lgd2,'Units', 'normalized','Position',[0.322 0.513972556141716 0.143 0.111270777479893],'FontSize',10.6);
set(lgd3,'Units', 'normalized','Position',[0.3055 0.258664245149759 0.159 0.0670241286863271],'FontSize',12)
set(lgd4,'Units', 'normalized','Position',[0.5705 0.710071751852172 0.1675 0.0670241286863271],'FontSize',12);
set(lgd5,'Units', 'normalized','Position',[0.76 0.540835826919197 0.145 0.0857908847184999],'FontSize',11);
set(lgd6,'Units', 'normalized','Position',[0.7375 0.260674969010349 0.1675 0.0650134048257373],'FontSize',12);

exportgraphics(f,'Fig5.png','Resolution', 400)
