
% Number of independent runs
S = 3;

% Number of parameters ranges from 2^2 to 2^l
l = 8; 
N = 2.^(3:l);
gram = false;
 
dataset = 'new';

list_cookies=[1,2,3,7,8,9];


gram =false;

if dataset == 'old'
if gram
    load("Cookies_gram_4_old.mat");
end
load("Cookies_randomized_then_orthogonalize_4_old.mat");
load("Cookies_randomized_then_orthogonalize_krp_4_old.mat");
load("Cookies_det_4_old.mat");
else
N = N(1:3);
if gram
    load("Cookies_nx100_gram_6_new.mat");
end
load("Cookies_nx100_randomized_then_orthogonalize_6_new.mat");
load("Cookies_nx100_randomized_then_orthogonalize_krp_6_new.mat");
load("Cookies_nx100_det_6_new.mat");
end
end_index = length(N);

% Discard results of first run
sumTimeNormalGMRES = sumTimeNormalGMRES(1:end_index,2:end);
otherTimeNormalGMRES = opTimeNormalGMRES+precTimeNormalGMRES+remTimeNormalGMRES;
otherTimeNormalGMRES = otherTimeNormalGMRES(1:end_index,2:end);
runtimeNormalGMRES = runtimeNormalGMRES(1:end_index,2:end);

if gram
    sumTimeGramGMRES = sumTimeGramGMRES(1:end_index,2:end);
    otherTimeGramGMRES = opTimeGramGMRES+precTimeGramGMRES+remTimeGramGMRES;
    otherTimeGramGMRES = otherTimeGramGMRES(1:end_index,2:end);
    runtimeGramGMRES = runtimeGramGMRES(1:end_index,2:end);
end

sumTimeRandGMRES = sumTimeRandGMRES(1:end_index,2:end);
otherTimeRandGMRES = opTimeRandGMRES+precTimeRandGMRES+remTimeRandGMRES;
otherTimeRandGMRES = otherTimeRandGMRES(1:end_index,2:end);
runtimeRandGMRES = runtimeRandGMRES(1:end_index,2:end);

sumTimeRandKRPGMRES = sumTimeRandKRPGMRES(1:end_index,2:end);
otherTimeRandKRPGMRES = opTimeRandKRPGMRES+precTimeRandKRPGMRES+remTimeRandKRPGMRES;
otherTimeRandKRPGMRES = otherTimeRandKRPGMRES(1:end_index,2:end);
runtimeRandKRPGMRES = runtimeRandKRPGMRES(1:end_index,2:end);

% if gram
%     save('test_case_with_gram_svd.mat', 'sumTimeNormalGMRES', 'otherTimeNormalGMRES', 'runtimeNormalGMRES','sumTimeGramGMRES', 'otherTimeGramGMRES', 'runtimeGramGMRES','sumTimeRandGMRES', 'otherTimeRandGMRES','runtimeRandGMRES','sumTimeRandKRPGMRES','otherTimeRandKRPGMRES','runtimeRandKRPGMRES','N','ranksrand','ranksnormal','ranksrandkrp',"ranksGram");
% else
%     save('test_case_without_gram_svd.mat', 'sumTimeNormalGMRES', 'otherTimeNormalGMRES', 'runtimeNormalGMRES','sumTimeRandGMRES', 'otherTimeRandGMRES','runtimeRandGMRES','sumTimeRandKRPGMRES','otherTimeRandKRPGMRES','runtimeRandKRPGMRES','N','ranksrand','ranksnormal','ranksrandkrp');
% end

%% Plot Results

close("all")

% load("test_case_without_gram_svd.mat")
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

if dataset == 'new'
file_name = "TTGMRES_Runtime_nx100_"+num2str(length(list_cookies))+"new.png";
else
file_name = "TTGMRES_Runtime_"+num2str(length(list_cookies))+"new.png";
end

if gram
    exportgraphics(f, 'TTGMRES_Runtime_with_Gram.png');
else
    exportgraphics(f, file_name)
end


f = figure(2);
f.Position(1:2) = [0,525];
f.Position(3:4) = [525, 350];

 
% %Deterministic with GRAM
% speedup = median(runtimeNormalGMRES, 2) ./ median(runtimeGramGMRES, 2);
% neg = speedup - min(runtimeNormalGMRES,[],2)./ max(runtimeGramGMRES,[],2);
% pos = max(runtimeNormalGMRES,[],2)./ min(runtimeGramGMRES,[],2) - speedup;
% errorbar(N, speedup, neg, pos,'c--s','markersize',10,'linewidth',2);
% hold on;
% speedup = median(sumTimeNormalGMRES, 2) ./ median(sumTimeGramGMRES, 2);
% neg = speedup - min(sumTimeNormalGMRES,[],2)./ max(sumTimeGramGMRES,[],2);
% pos = max(sumTimeNormalGMRES,[],2)./ min(sumTimeGramGMRES,[],2) - speedup;
% errorbar(N, speedup, neg, pos,'k--o','markersize',10,'linewidth',2);

if gram
    %Randomized with TT Structure
    speedup = median(runtimeGramGMRES, 2) ./ median(runtimeRandGMRES, 2);
    neg = speedup - min(runtimeGramGMRES,[],2)./ max(runtimeRandGMRES,[],2);
    pos = max(runtimeGramGMRES,[],2)./ min(runtimeRandGMRES,[],2) - speedup;
    errorbar(N, speedup, neg, pos,'r--s','markersize',10,'linewidth',2);
    hold on;
    speedup = median(sumTimeGramGMRES, 2) ./ median(sumTimeRandGMRES, 2);
    neg = speedup - min(sumTimeGramGMRES,[],2)./ max(sumTimeRandGMRES,[],2);
    pos = max(sumTimeGramGMRES,[],2)./ min(sumTimeRandGMRES,[],2) - speedup;
    errorbar(N, speedup, neg, pos,'b--o','markersize',10,'linewidth',2);
end


%Randomized with TT Structure
speedup = median(runtimeNormalGMRES, 2) ./ median(runtimeRandGMRES, 2);
neg = speedup - min(runtimeNormalGMRES,[],2)./ max(runtimeRandGMRES,[],2);
pos = max(runtimeNormalGMRES,[],2)./ min(runtimeRandGMRES,[],2) - speedup;
errorbar(N, speedup, neg, pos,'-s','Color','#7E2F8E','markersize',10,'linewidth',2);
hold on;
speedup = median(sumTimeNormalGMRES, 2) ./ median(sumTimeRandGMRES, 2);
neg = speedup - min(sumTimeNormalGMRES,[],2)./ max(sumTimeRandGMRES,[],2);
pos = max(sumTimeNormalGMRES,[],2)./ min(sumTimeRandGMRES,[],2) - speedup;
errorbar(N, speedup, neg, pos,'-o','Color','#008000','markersize',10,'linewidth',2);

if gram
    % %Randomized with KRP Structure
    % speedup_krp = median(runtimeGramGMRES, 2) ./ median(runtimeRandKRPGMRES, 2);
    % neg_krp = speedup_krp - min(runtimeGramGMRES,[],2)./ max(runtimeRandKRPGMRES,[],2);
    % pos_krp = max(runtimeGramGMRES,[],2)./ min(runtimeRandKRPGMRES,[],2) - speedup_krp;
    % errorbar(N, speedup_krp, neg_krp, pos_krp,'m--s','markersize',10,'linewidth',2);
    % speedup_krp = median(sumTimeGramGMRES, 2) ./ median(sumTimeRandKRPGMRES, 2);
    % neg_krp = speedup_krp - min(sumTimeGramGMRES,[],2)./ max(sumTimeRandKRPGMRES,[],2);
    % pos_krp = max(sumTimeGramGMRES,[],2)./ min(sumTimeRandKRPGMRES,[],2) - speedup_krp;
    % errorbar(N, speedup_krp, neg_krp, pos_krp,'g--o','markersize',10,'linewidth',2);
end

%Randomized with KRP Structure
speedup_krp = median(runtimeNormalGMRES, 2) ./ median(runtimeRandKRPGMRES, 2);
neg_krp = speedup_krp - min(runtimeNormalGMRES,[],2)./ max(runtimeRandKRPGMRES,[],2);
pos_krp = max(runtimeNormalGMRES,[],2)./ min(runtimeRandKRPGMRES,[],2) - speedup_krp;
errorbar(N, speedup_krp, neg_krp, pos_krp,'r-s','markersize',10,'linewidth',2);
speedup_krp = median(sumTimeNormalGMRES, 2) ./ median(sumTimeRandKRPGMRES, 2);
neg_krp = speedup_krp - min(sumTimeNormalGMRES,[],2)./ max(sumTimeRandKRPGMRES,[],2);
pos_krp = max(sumTimeNormalGMRES,[],2)./ min(sumTimeRandKRPGMRES,[],2) - speedup_krp;
errorbar(N, speedup_krp, neg_krp, pos_krp,'b-o','markersize',10,'linewidth',2);
hold off;

if gram
    legend('Total speedup-Gram', 'TT-Sum + Round-Gram','Total speedup', 'TT-Sum + Round','Total speedup--KRP','TT-Sum + Round--KRP', 'Total speedup-KRP','TT-Sum + Round-KRP', 'location', 'northwest')
else
    legend('Total speedup-TTlike', 'TT-Sum+Round-TTlike','Total speedup-KRP','TT-Sum+Round-KRP', 'location', 'southeast')
end

set(gca, 'XScale', 'log')
xlabel('Number of parameter samples I_2 = \ldots = I_N (log scale)', 'FontSize', 18)
ylabel('Speedup', 'FontSize', 18)
xticks(N);
set(gca,'FontSize',16)
% grid on

if dataset == 'new'
    file_name = "TTGMRES_Speedup_nx100_"+num2str(length(list_cookies))+"new.png";
else
    file_name = "TTGMRES_Speedup_"+num2str(length(list_cookies))+"new.png";
end
if gram
    exportgraphics(f, file_name)
else
    exportgraphics(f, file_name)
    % exportgraphics(f, 'TTGMRES_Speedup.pdf')
end

f = figure(3);
clf();
f.Position(1:2) = [525,525];
f.Position(3:4) = [525, 355];

ranks_krp = max(ranksrandkrp{1}(1:end-1,:),[],2);
p3 = plot(ranks_krp, 'r*','markersize',10,'linewidth',2); hold on

ranks_rand = max(ranksrand{1}(1:end-1,:),[],2);
p1 = plot(ranks_rand, 'bo','markersize',10,'linewidth',2);

ranks_det = max(ranksnormal{1}(1:end-1,:),[],2);
p2 = plot(ranks_det, 'k+','markersize',10,'linewidth',2);


if gram
    ranks = max(ranksGram{1}(1:end-1,:),[],2);
    p4 = plot(ranks, 'm*','markersize',10,'linewidth',2);
end

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

if gram
    for j=2:length(N)
        ranks = max(ranksGram{j}(1:end-1,:),[],2);
        plot(ranks, 'm*', 'markersize',10,'linewidth',2)
    end
end


hold off
if gram
    legend([p2,p4,p1,p3], 'Naive TT-GMRES', 'Gram TT-GMRES','Randomized TT-GMRES-TTlike','Randomized TT-GMRES-KRP' , 'location', 'southeast')
else
    legend([p3,p1,p2], 'Randomized TT-GMRES-KRP' , 'Randomized TT-GMRES-TTlike','Naive TT-GMRES', 'location', 'southeast')

    % legend([p2,p1,p3], 'Naive TT-GMRES', 'Randomized TT-GMRES-TTlike','Randomized TT-GMRES-KRP' , 'location', 'southeast')
end
xlabel('TT-GMRES iteration number', 'FontSize', 18)
ylabel('TT-rank', 'FontSize', 18)
set(gca,'FontSize',16)

if dataset == 'new'
    file_name = "TTGMRES_ranks_nx100_"+num2str(length(list_cookies))+"new.png";
else
    file_name = "TTGMRES_ranks_"+num2str(length(list_cookies))+"new.png";
end

if gram
    exportgraphics(f, file_name)
else
    exportgraphics(f, file_name)
end

beep;
pause(1)
beep;