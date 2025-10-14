clc
clear all 
close all

load('Synthetic_results.mat')

%% Plot Results

% Post-process errors
[errLR,       negLR,       posLR]       = computeError(errorsLR);
[errRandLR,   negRandLR,   posRandLR]   = computeError(errorsRandLR);
[errRandOrth, negRandOrth, posRandOrth] = computeError(errorsRandOrth);
[errRandOrth_krp, negRandOrth_krp, posRandOrth_krp] = computeError(errorsRandOrth_krp);

f = figure('Position', [100, 100, 800, 280]);
subplot(1,2,1)
hold on
p1 = errorbar(test_ranks, errLR, negLR, posLR,     'ko','markersize',10,'linewidth',2,'DisplayName','TT-Rounding');
p2 = errorbar(test_ranks, errRandLR,   negRandLR,   posRandLR,    'b+','markersize',10,'linewidth',2,'DisplayName','TT-Orth-Rand (Fix)');
p4 = errorbar(test_ranks, errRandOrth_krp, negRandOrth_krp, posRandOrth_krp,  'r*','markersize',10,'linewidth',2, 'DisplayName','TT-KRP (Fix)');
p3 = errorbar(test_ranks, errRandOrth, negRandOrth, posRandOrth,  'x','Color','#7E2F8E','markersize',10,'linewidth',2,'DisplayName', 'TT-Rand-Orth');
hold off

yticks([1e-5,1e-4,1e-3,1e-2,1e-1, 1e0])
ylim([5e-6,2e0])
ax = gca;
ax.YScale = 'log';
xlabel('Maximum Target Rank', 'FontSize', 16)
ylabel('Relative Error', 'FontSize', 16)
legend([p4,p3,p2,p1],'Location', 'northeast')
set(gca,'FontSize',15)
grid on;
box on;


% Post-process timings

[speedupLR,         negLR,          posLR]          = computeSpeedup(timeLR,       timeLR);
[speedupRandLR,     negRandLR,      posRandLR]      = computeSpeedup(timeRandLR,   timeLR);
[speedupRandOrth,   negRandOrth,    posRandOrth]    = computeSpeedup(timeRandOrth, timeLR);
[speedupRandOrth_krp,   negRandOrth_krp,    posRandOrth_krp]    = computeSpeedup(timeRandOrth_krp, timeLR);

subplot(1,2,2)
hold on
p1= plot(test_ranks, speedupLR,       'ko-','markersize',10,'linewidth',2,'DisplayName','TT-Rounding');
p2 = plot(test_ranks, speedupRandLR,     'b+-','markersize',10,'linewidth',2,'DisplayName','TT-Orth-Rand (Fix)');
p3 = plot(test_ranks, speedupRandOrth,  'x-','Color','#7E2F8E','markersize',10,'linewidth',2,'DisplayName', 'TT-Rand-Orth');
p4 = plot(test_ranks, speedupRandOrth_krp, 'r*-','markersize',10,'linewidth',2,'DisplayName','TT-KRP (Fix)');
    
hold off
xlabel('Maximum Target Rank', 'FontSize', 16)
ylabel('Speedup', 'FontSize', 16)
legend([p4,p3,p2,p1],'Location', 'northeast')
set(gca,'FontSize',15)
box on;
grid on;

exportgraphics(f, 'Fig4.png','Resolution', 600);



%% Custom function to post-process errors

function [error, neg, pos] = computeError(errors)
    errors = errors(:,2:end);
    error = squeeze(median(errors, 2));
    neg = error -  squeeze(min(errors,[],2));
    pos = squeeze(max(errors,[],2)) - error;
end

%% Custom function to post-process runtimes

function [speedup, neg, pos] = computeSpeedup(time, timeref)
    time = time(:,2:end);
    timeref = timeref(:,2:end);
    speedup = mean(timeref,2) ./ mean(time,2);
    neg = speedup - median(timeref,2) ./ max(time,[],2);
    pos = median(timeref,2) ./ min(time,[],2) - speedup;
end