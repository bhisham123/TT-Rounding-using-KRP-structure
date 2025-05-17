%% TTGMRES test
clear;
clc;
maxNumCompThreads(2);
close("all")


%% Add path
addpath("./TTcore/");
addpath("./TTrandomized/");
addpath("./Data/");
addpath("./Cookies/");
addpath("./Gram/")


%% Global parameters

% Number of independent runs
S = 3;

% Number of parameters ranges from 2^2 to 2^l
l = 6; %8; 
N = 2.^(3:l);


%% Load FEM matrices
gram = false;
A0 = readcoomat('A0.txt');
A1 = readcoomat('A1.txt');
A2 = readcoomat('A2.txt');
A3 = readcoomat('A3.txt');
A4 = readcoomat('A4.txt');

n1 = length(A0);
A = speye(n1);

a0 = readb('b.txt');

%% Run and time experiments

sumTimeNormalGMRES = zeros(length(N),S);
sumTimeGramGMRES = zeros(length(N),S);
sumTimeRandGMRES = zeros(length(N),S);
sumTimeRandKRPGMRES = zeros(length(N),S);


opTimeNormalGMRES = zeros(length(N),S);
opTimeGramGMRES = zeros(length(N),S);
opTimeRandGMRES = zeros(length(N),S);
opTimeRandKRPGMRES = zeros(length(N),S);

precTimeNormalGMRES = zeros(length(N),S);
precTimeGramGMRES = zeros(length(N),S);
precTimeRandGMRES = zeros(length(N),S);
precTimeRandKRPGMRES = zeros(length(N),S);


remTimeNormalGMRES = zeros(length(N),S);
remTimeGramGMRES = zeros(length(N),S);
remTimeRandGMRES = zeros(length(N),S);
remTimeRandKRPGMRES = zeros(length(N),S);


runtimeNormalGMRES = zeros(length(N),S);
runtimeGramGMRES = zeros(length(N),S);
runtimeRandGMRES = zeros(length(N),S);
runtimeRandKRPGMRES = zeros(length(N),S);

ranksrand = cell(length(N));
ranksrandkrp = cell(length(N));
ranksnormal = cell(length(N));
ranksGram = cell(length(N));


for i = 1:length(N)
    n = N(i);

    % Distribution of parameters
    MinD = 1e0;
    MaxD = 1e1;

    d = linspace(MinD, MaxD, n);

    % Matrices and operator in Kronecker form
    D = sparse(1:n, 1:n, d);

    c = ones(n,1);
    B = speye(n);

    b = {a0;c;c;c;c};

    A = {A0, A1, A2, A3, A4;
         B,  D,  B,  B,  B;
         B,  B,  D,  B,  B;
         B,  B,  B,  D,  B;
         B,  B,  B,  B,  D};

    % TT-GMRES parameter setup
    tol      = 1e-8;
    maxit    = 50;
    [L,U,p]  = lu(A0 + 1*(A1 + A2 + A3 + A4), 'vector');
    prec     = @(x) Preconditioner(L,U,p,x);
    tols     = tol*1e-2;
    normb    = TTnorm(b);
    b        = TTscale(b, 1/normb);
    normb    = TTnorm(b);

    Op            = @(x)       TTsummandsKronOp(A, x);
    
    deterministic = @(W,a,tol) TTsum(W, a, tol);
    deterministic_gram = @(W,a,tol) TTsum_Gram(W, a, tol);
    randomized    = @(W,a,tol) TTsum_Randomize_then_Orthogonalize(W, a, tol);
    randomized_krp    = @(W,a,tol) TTsum_Randomize_then_Orthogonalize_KRP(W, a, tol);

    for s=1:S

        if gram
            fprintf("*********\nDeterministic GRAM TT-GMRES run - n2 = %i\t", n);
            tStart = tic;
            [x, ranksN, sumtime, optime, prectime, remtime] = timed_TTGMRES(Op, b, tol, maxit, prec, tols, deterministic_gram);
            tEnd = toc(tStart);
            fprintf("Elapsed time is %.2f seconds.\n", tEnd);
            t = TTsummandsKronOp(A,x);
            r = t{1};
            for j=2:length(t)
                r = TTaxby(1., r, 1., t{j});
            end
            r = TTaxby(1., r, -1., b);
            fprintf("*********\nnorm real res = %e\t", TTnorm(r, "OLR")/normb);
            fprintf("rounding and sum time = %.2f\n*********\n", sumtime);
    
            for j = 1 : size(ranksN,1)
                ranksGram{i}(j,s) = max(ranksN{j}); 
            end
            sumTimeGramGMRES(i,s)  = sumtime;
            opTimeGramGMRES(i,s)   = optime;
            precTimeGramGMRES(i,s) = prectime;
            remTimeGramGMRES(i,s)  = remtime;
            runtimeGramGMRES(i,s)  = tEnd;
        end

    
        fprintf("*********\nRandomized TT-GMRES run - n2 = %i\t", n);
        tStart = tic;
        [x, ranksR, sumtime, optime, prectime, remtime] = timed_TTGMRES(Op, b, tol, maxit, prec, tols, randomized);
        tEnd = toc(tStart);
        fprintf("Elapsed time is %.2f seconds.\n", tEnd);
        t = TTsummandsKronOp(A,x);
        r = t{1};
        for j=2:length(t)
            r = TTaxby(1., r, 1., t{j});
        end
        r = TTaxby(1., r, -1., b);
        fprintf("*********\nnorm real res = %e\t", TTnorm(r, "OLR")/normb);
        fprintf("rounding and sum time = %.2f\n********\n", sumtime);

        for j = 1 : size(ranksR,1)
            ranksrand{i}(j,s) = max(ranksR{j});
        end
        sumTimeRandGMRES(i,s)   = sumtime;
        opTimeRandGMRES(i,s)    = optime;
        precTimeRandGMRES(i,s)  = prectime;
        remTimeRandGMRES(i,s)   = remtime;
        runtimeRandGMRES(i,s)   = tEnd;


        fprintf("*********\nRandomized-KRP TT-GMRES run - n2 = %i\t", n);
        tStart = tic;
        [x, ranksR_KRP, sumtime, optime, prectime, remtime] = timed_TTGMRES(Op, b, tol, maxit, prec, tols, randomized_krp);
        tEnd = toc(tStart);
        fprintf("Elapsed time is %.2f seconds.\n", tEnd);
        t = TTsummandsKronOp(A,x);
        r = t{1};
        for j=2:length(t)
            r = TTaxby(1., r, 1., t{j});
        end
        r = TTaxby(1., r, -1., b);
        fprintf("*********\nnorm real res = %e\t", TTnorm(r, "OLR")/normb);
        fprintf("rounding and sum time = %.2f\n********\n", sumtime);

        for j = 1 : size(ranksR_KRP,1)
            ranksrandkrp{i}(j,s) = max(ranksR_KRP{j});
        end
        sumTimeRandKRPGMRES(i,s)   = sumtime;
        opTimeRandKRPGMRES(i,s)    = optime;
        precTimeRandKRPGMRES(i,s)  = prectime;
        remTimeRandKRPGMRES(i,s)   = remtime;
        runtimeRandKRPGMRES(i,s)   = tEnd;

        
        fprintf("*********\nDeterministic TT-GMRES run - n2 = %i\t", n);
        tStart = tic;
        [x, ranksN, sumtime, optime, prectime, remtime] = timed_TTGMRES(Op, b, tol, maxit, prec, tols, deterministic);
        tEnd = toc(tStart);
        fprintf("Elapsed time is %.2f seconds.\n", tEnd);
        t = TTsummandsKronOp(A,x);
        r = t{1};
        for j=2:length(t)
            r = TTaxby(1., r, 1., t{j});
        end
        r = TTaxby(1., r, -1., b);
        fprintf("*********\nnorm real res = %e\t", TTnorm(r, "OLR")/normb);
        fprintf("rounding and sum time = %.2f\n*********\n", sumtime);
        
        for j = 1 : size(ranksN,1)
            ranksnormal{i}(j,s) = max(ranksN{j});
        end
        sumTimeNormalGMRES(i,s)  = sumtime;
        opTimeNormalGMRES(i,s)   = optime;
        precTimeNormalGMRES(i,s) = prectime;
        remTimeNormalGMRES(i,s)  = remtime;
        runtimeNormalGMRES(i,s)  = tEnd;   
    end
end

% Discard results of first run
sumTimeNormalGMRES = sumTimeNormalGMRES(:,2:end);
otherTimeNormalGMRES = opTimeNormalGMRES+precTimeNormalGMRES+remTimeNormalGMRES;
otherTimeNormalGMRES = otherTimeNormalGMRES(:,2:end);
runtimeNormalGMRES = runtimeNormalGMRES(:,2:end);

sumTimeGramGMRES = sumTimeGramGMRES(:,2:end);
otherTimeGramGMRES = opTimeGramGMRES+precTimeGramGMRES+remTimeGramGMRES;
otherTimeGramGMRES = otherTimeGramGMRES(:,2:end);
runtimeGramGMRES = runtimeGramGMRES(:,2:end);


sumTimeRandGMRES = sumTimeRandGMRES(:,2:end);
otherTimeRandGMRES = opTimeRandGMRES+precTimeRandGMRES+remTimeRandGMRES;
otherTimeRandGMRES = otherTimeRandGMRES(:,2:end);
runtimeRandGMRES = runtimeRandGMRES(:,2:end);

sumTimeRandKRPGMRES = sumTimeRandKRPGMRES(:,2:end);
otherTimeRandKRPGMRES = opTimeRandKRPGMRES+precTimeRandKRPGMRES+remTimeRandKRPGMRES;
otherTimeRandKRPGMRES = otherTimeRandKRPGMRES(:,2:end);
runtimeRandKRPGMRES = runtimeRandKRPGMRES(:,2:end);

% if gram
%     save('test_case_with_gram_svd.mat', 'sumTimeNormalGMRES', 'otherTimeNormalGMRES', 'runtimeNormalGMRES','sumTimeGramGMRES', 'otherTimeGramGMRES', 'runtimeGramGMRES','sumTimeRandGMRES', 'otherTimeRandGMRES','runtimeRandGMRES','sumTimeRandKRPGMRES','otherTimeRandKRPGMRES','runtimeRandKRPGMRES','N','ranksrand','ranksnormal','ranksrandkrp',"ranksGram");
% else
%     save('test_case_without_gram_svd.mat', 'sumTimeNormalGMRES', 'otherTimeNormalGMRES', 'runtimeNormalGMRES','sumTimeRandGMRES', 'otherTimeRandGMRES','runtimeRandGMRES','sumTimeRandKRPGMRES','otherTimeRandKRPGMRES','runtimeRandKRPGMRES','N','ranksrand','ranksnormal','ranksrandkrp');
% end

%% Plot Results

close("all")

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

% if gram
%     exportgraphics(f, 'TTGMRES_Runtime_with_Gram.png')
% else
%     exportgraphics(f, 'TTGMRES_Runtime.png')
% end


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
errorbar(N, speedup, neg, pos,'r-s','markersize',10,'linewidth',2);
hold on;
speedup = median(sumTimeNormalGMRES, 2) ./ median(sumTimeRandGMRES, 2);
neg = speedup - min(sumTimeNormalGMRES,[],2)./ max(sumTimeRandGMRES,[],2);
pos = max(sumTimeNormalGMRES,[],2)./ min(sumTimeRandGMRES,[],2) - speedup;
errorbar(N, speedup, neg, pos,'b-o','markersize',10,'linewidth',2);

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
errorbar(N, speedup_krp, neg_krp, pos_krp,'m-s','markersize',10,'linewidth',2);
speedup_krp = median(sumTimeNormalGMRES, 2) ./ median(sumTimeRandKRPGMRES, 2);
neg_krp = speedup_krp - min(sumTimeNormalGMRES,[],2)./ max(sumTimeRandKRPGMRES,[],2);
pos_krp = max(sumTimeNormalGMRES,[],2)./ min(sumTimeRandKRPGMRES,[],2) - speedup_krp;
errorbar(N, speedup_krp, neg_krp, pos_krp,'g-o','markersize',10,'linewidth',2);
hold off;

if gram
    legend('Total speedup-Gram', 'TT-Sum + Round-Gram','Total speedup', 'TT-Sum + Round','Total speedup--KRP','TT-Sum + Round--KRP', 'Total speedup-KRP','TT-Sum + Round-KRP', 'location', 'northwest')
else
    legend('Total speedup', 'TT-Sum + Round','Total speedup--KRP','TT-Sum + Round--KRP', 'location', 'northwest')
end

set(gca, 'XScale', 'log')
xlabel('Number of parameter samples I_2 = \ldots = I_N (log scale)', 'FontSize', 18)
ylabel('Speedup', 'FontSize', 18)
xticks(N);
set(gca,'FontSize',16)
grid on

% if gram
%     exportgraphics(f, 'TTGMRES_Speedup_with_Gram.png')
% else
%     exportgraphics(f, 'TTGMRES_Speedup.png')
% end

f = figure(3);
clf();
f.Position(1:2) = [525,525];
f.Position(3:4) = [525, 355];
ranks = max(ranksrand{1}(1:end-1,:),[],2);
p1 = plot(ranks, 'bo','markersize',10,'linewidth',2);
hold on
ranks = max(ranksnormal{1}(1:end-1,:),[],2);
p2 = plot(ranks, 'r+','markersize',10,'linewidth',2);

ranks = max(ranksrandkrp{1}(1:end-1,:),[],2);
p3 = plot(ranks, 'ks','markersize',10,'linewidth',2);

if gram
    ranks = max(ranksGram{1}(1:end-1,:),[],2);
    p4 = plot(ranks, 'm*','markersize',10,'linewidth',2);
end

for j=2:length(N)
    ranks = max(ranksrand{j}(1:end-1,:),[],2);
    plot(ranks, 'bo', 'markersize',10,'linewidth',2)
end
for j=2:length(N)
    ranks = max(ranksnormal{j}(1:end-1,:),[],2);
    plot(ranks, 'r+', 'markersize',10,'linewidth',2)
end

if gram
    for j=2:length(N)
        ranks = max(ranksGram{j}(1:end-1,:),[],2);
        plot(ranks, 'm*', 'markersize',10,'linewidth',2)
    end
end

for j=2:length(N)
    ranks = max(ranksrandkrp{j}(1:end-1,:),[],2);
    plot(ranks, 'ks', 'markersize',10,'linewidth',2)
end
hold off
if gram
    legend([p2,p4,p1,p3], 'Naive TT-GMRES', 'Gram TT-GMRES','Randomized TT-GMRES','Randomized TT-GMRES-KRP' , 'location', 'southeast')
else
    legend([p2,p1,p3], 'Naive TT-GMRES', 'Randomized TT-GMRES','Randomized TT-GMRES-KRP' , 'location', 'southeast')
end
xlabel('TT-GMRES iteration number', 'FontSize', 18)
ylabel('TT-rank', 'FontSize', 18)
set(gca,'FontSize',16)

% if gram
%     exportgraphics(f, 'TTGMRES_Ranks_with_Gram.png')
% else
%     exportgraphics(f, 'TTGMRES_Ranks.png')
% end

beep;
pause(1)
beep;


%% Preconditioner
function X = Preconditioner(L,U,p, Y)
    X = Y;
    X{1} = U\(L\Y{1}(p,:));
end
