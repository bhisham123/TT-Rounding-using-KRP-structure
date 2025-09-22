%% TTGMRES test
clear;
clc;
maxNumCompThreads(2);
close("all")


%% Add path
addpath("./TTcore/");
addpath("./TTrandomized/");
addpath("./Data/");
addpath("./cookies/");
addpath("./Gram/");
nx = 'nx_40';
addpath("./Data/nine_cookies/"+nx+"/");


%% Global parameters

% Number of independent runs
S = 3;

% Number of parameters ranges from 2^2 to 2^l
l = 8; 
N = 2.^(3:l);
gram = false;

% select the cookies data 
dataset = 'new';

% select a number of the nine cookies
list_cookies=[1,2,3,5,7,8,9];

%% Load FEM matrices
if dataset == 'old'
    list_cookies=[1,2,3,4];
    A{1} = readcoomat('A0.txt');
    A{2} = readcoomat('A1.txt');
    A{3} = readcoomat('A2.txt');
    A{4} = readcoomat('A3.txt');
    A{5} = readcoomat('A4.txt');    
    a0 = readb('b.txt');
    A0 = A{1};
else dataset == 'new'
    matrix_1;
    matrix_2;
    matrix_3;
    matrix_4;
    matrix_5;
    matrix_6;
    matrix_7;
    matrix_8;
    matrix_9;
    matrix_10;
    rhs_;
    a0 = rhs;
    clear rhs;
    A0 = A{1};
    size(A0)
end

n1 = length(A0);
% A = speye(n1);



%% Run and time experiments

sumTimeNormalGMRES = zeros(length(N),S);
sumTimeGramGMRES = zeros(length(N),S);
sumTimeRandGMRES = zeros(length(N),S);
sumTimeRandKRPGMRES = zeros(length(N),S);

res_err_normal = zeros(length(N),S);
res_err_gram = zeros(length(N),S);
res_err_rand = zeros(length(N),S);
res_err_randKRP = zeros(length(N),S);


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
    MaxD = 5; %1e1;

    d = linspace(MinD, MaxD, n);

    % Matrices and operator in Kronecker form
    D = sparse(1:n, 1:n, d);

    c = ones(n,1);
    B = speye(n);
    b = cell(length(list_cookies)+1,1);
    b{1} = a0;
    for j = 1 : length(list_cookies)
        b{j+1} = c;
    end
    C = cell(length(list_cookies)+1,length(list_cookies)+1);
    for p = 1 : length(list_cookies) + 1
        for q = 1 : length(list_cookies) + 1
            if p == 1
                % C{1,q} = A{p};
                if q == 1
                    C{1,q} = A{q};
                    fprintf('A{%d}\t',q);
                else
                    C{1,q} = A{list_cookies(q-1)+1};
                    fprintf('A{%d}\t',list_cookies(q-1)+1);
                end

            else
                if p == q
                    C{p,q} = D;
                    fprintf('D{0}\t');
                else
                    C{p,q} = B;
                    fprintf('B{0}\t')
                end
            end
        end
        fprintf('\n');
    end
    % A = C;
    % clear C;
    

    % TT-GMRES parameter setup
    tol      = 1e-8;
    maxit    = 50;

    mat = C{1,1};
    for j = 1 : length(list_cookies)
        % mat = mat + A{1,i+1};
        mat = mat + C{1,j+1};
        fprintf('A{1,%d} + ',list_cookies(j)+1);
    end
    [L,U,p]  = lu(mat, 'vector');
    clear mat;
    prec     = @(x) Preconditioner(L,U,p,x);
    tols     = tol*1e-2;
    normb    = TTnorm(b);
    b        = TTscale(b, 1/normb);
    normb    = TTnorm(b);

    Op            = @(x)       TTsummandsKronOp(C, x);

    deterministic = @(W,a,tol) TTsum(W, a, tol);
    deterministic_gram = @(W,a,tol) TTsum_Gram(W, a, tol);
    randomized    = @(W,a,tol) TTsum_Randomize_then_Orthogonalize(W, a, tol);
    randomized_krp    = @(W,a,tol) TTsum_Randomize_then_Orthogonalize_KRP(W, a, tol);

    if gram
        for s=1:S
            fprintf("\n*********\nDeterministic GRAM TT-GMRES run - n2 = %i\t", n);
            tStart = tic;
            [x, ranksN, sumtime, optime, prectime, remtime] = timed_TTGMRES(Op, b, tol, maxit, prec, tols, deterministic_gram);
            tEnd = toc(tStart);
            fprintf("Elapsed time is %.2f seconds.\n", tEnd);
            t = TTsummandsKronOp(C,x);
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

            res_err_gram(i,s) = TTnorm(r, "OLR")/normb;
            sumTimeGramGMRES(i,s)  = sumtime;
            opTimeGramGMRES(i,s)   = optime;
            precTimeGramGMRES(i,s) = prectime;
            remTimeGramGMRES(i,s)  = remtime;
            runtimeGramGMRES(i,s)  = tEnd;

            file_name = "Cookies_gram_"+num2str(length(list_cookies))+"_"+dataset+"_"+nx+"_maxD_"+num2str(MaxD)+".mat";
            % save(file_name, 'sumTimeGramGMRES', 'opTimeGramGMRES','precTimeGramGMRES','remTimeGramGMRES','runtimeGramGMRES','ranksGram','gram','N','dataset','list_cookies','res_err_gram');
        end
    end

    for s=1:S
        fprintf("\n*********\nRandomized TT-GMRES run - n2 = %i\t", n);
        tStart = tic;
        [x, ranksR, sumtime, optime, prectime, remtime] = timed_TTGMRES(Op, b, tol, maxit, prec, tols, randomized);
        tEnd = toc(tStart);
        fprintf("Elapsed time is %.2f seconds.\n", tEnd);
        t = TTsummandsKronOp(C,x);
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

        res_err_rand(i,s) = TTnorm(r, "OLR")/normb;
        sumTimeRandGMRES(i,s)   = sumtime;
        opTimeRandGMRES(i,s)    = optime;
        precTimeRandGMRES(i,s)  = prectime;
        remTimeRandGMRES(i,s)   = remtime;
        runtimeRandGMRES(i,s)   = tEnd;

        file_name = "Cookies_randomized_then_orthogonalize_"+num2str(length(list_cookies))+"_"+dataset+"_"+nx+"_maxD_"+num2str(MaxD)+".mat";
        % save(file_name, 'sumTimeRandGMRES', 'opTimeRandGMRES','precTimeRandGMRES','remTimeRandGMRES','runtimeRandGMRES','ranksrand','gram','N','dataset','list_cookies','res_err_rand');
    end

    for s=1:S
        fprintf("*********\nRandomized-KRP TT-GMRES run - n2 = %i\t", n);
        tStart = tic;
        [x, ranksR_KRP, sumtime, optime, prectime, remtime] = timed_TTGMRES(Op, b, tol, maxit, prec, tols, randomized_krp);
        tEnd = toc(tStart);
        fprintf("Elapsed time is %.2f seconds.\n", tEnd);
        t = TTsummandsKronOp(C,x);
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

        res_err_randKRP(i,s) = TTnorm(r, "OLR")/normb;
        sumTimeRandKRPGMRES(i,s)   = sumtime;
        opTimeRandKRPGMRES(i,s)    = optime;
        precTimeRandKRPGMRES(i,s)  = prectime;
        remTimeRandKRPGMRES(i,s)   = remtime;
        runtimeRandKRPGMRES(i,s)   = tEnd;

        file_name = "Cookies_randomized_then_orthogonalize_krp_"+num2str(length(list_cookies))+"_"+dataset+"_"+nx+"_maxD_"+num2str(MaxD)+".mat";
        % save(file_name, 'sumTimeRandKRPGMRES', 'opTimeRandKRPGMRES','precTimeRandKRPGMRES','remTimeRandKRPGMRES','runtimeRandKRPGMRES','ranksrandkrp','gram','N','dataset','list_cookies','res_err_randKRP');
    end

    for s=1:S
        fprintf("*********\nDeterministic TT-GMRES run - n2 = %i\t", n);
        tStart = tic;
        [x, ranksN, sumtime, optime, prectime, remtime] = timed_TTGMRES(Op, b, tol, maxit, prec, tols, deterministic);
        tEnd = toc(tStart);
        fprintf("Elapsed time is %.2f seconds.\n", tEnd);
        t = TTsummandsKronOp(C,x);
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

        res_err_normal(i,s) = TTnorm(r, "OLR")/normb;
        sumTimeNormalGMRES(i,s)  = sumtime;
        opTimeNormalGMRES(i,s)   = optime;
        precTimeNormalGMRES(i,s) = prectime;
        remTimeNormalGMRES(i,s)  = remtime;
        runtimeNormalGMRES(i,s)  = tEnd;   


        file_name = "Cookies_det_"+num2str(length(list_cookies))+"_"+dataset+"_"+nx+"_maxD_"+num2str(MaxD)+".mat";
        % save(file_name, 'sumTimeNormalGMRES', 'opTimeNormalGMRES','precTimeNormalGMRES','remTimeNormalGMRES','runtimeNormalGMRES','ranksnormal','gram','N','dataset','list_cookies','res_err_normal');
    end
end
% clear all
% % Number of parameters ranges from 2^2 to 2^l
% l = 4; %8; 
% N = 2.^(3:l);
% gram = false;
% 
% dataset = '9';
% list_cookies=[1,3,7,9];

% if gram
%     load("Cookies_gram_"+num2str(length(list_cookies))+"_"+dataset+"_"+nx+".mat");
% end
% load("Cookies_randomized_then_orthogonalize_"+num2str(length(list_cookies))+"_"+dataset+"_"+nx+".mat");
% load("Cookies_randomized_then_orthogonalize_krp_"+num2str(length(list_cookies))+"_"+dataset+"_"+nx+".mat");
% load("Cookies_det_"+num2str(length(list_cookies))+"_"+dataset+"_"+nx+".mat");

% Discard results of first run
sumTimeNormalGMRES = sumTimeNormalGMRES(:,2:end);
otherTimeNormalGMRES = opTimeNormalGMRES+precTimeNormalGMRES+remTimeNormalGMRES;
otherTimeNormalGMRES = otherTimeNormalGMRES(:,2:end);
runtimeNormalGMRES = runtimeNormalGMRES(:,2:end);

if gram
    sumTimeGramGMRES = sumTimeGramGMRES(:,2:end);
    otherTimeGramGMRES = opTimeGramGMRES+precTimeGramGMRES+remTimeGramGMRES;
    otherTimeGramGMRES = otherTimeGramGMRES(:,2:end);
    runtimeGramGMRES = runtimeGramGMRES(:,2:end);
end

sumTimeRandGMRES = sumTimeRandGMRES(:,2:end);
otherTimeRandGMRES = opTimeRandGMRES+precTimeRandGMRES+remTimeRandGMRES;
otherTimeRandGMRES = otherTimeRandGMRES(:,2:end);
runtimeRandGMRES = runtimeRandGMRES(:,2:end);

sumTimeRandKRPGMRES = sumTimeRandKRPGMRES(:,2:end);
otherTimeRandKRPGMRES = opTimeRandKRPGMRES+precTimeRandKRPGMRES+remTimeRandKRPGMRES;
otherTimeRandKRPGMRES = otherTimeRandKRPGMRES(:,2:end);
runtimeRandKRPGMRES = runtimeRandKRPGMRES(:,2:end);

% file_name = "Cookies_"+num2str(length(list_cookies))+"_"+dataset+"_"+nx+"_maxD_"+num2str(MaxD)+".mat";
% if gram
%     save('test_case_with_gram_svd.mat', 'sumTimeNormalGMRES', 'otherTimeNormalGMRES', 'runtimeNormalGMRES','sumTimeGramGMRES', 'otherTimeGramGMRES', 'runtimeGramGMRES','sumTimeRandGMRES', 'otherTimeRandGMRES','runtimeRandGMRES','sumTimeRandKRPGMRES','otherTimeRandKRPGMRES','runtimeRandKRPGMRES','N','ranksrand','ranksnormal','ranksrandkrp',"ranksGram");
% else
%     save(file_name, 'sumTimeNormalGMRES', 'otherTimeNormalGMRES', 'runtimeNormalGMRES','sumTimeRandGMRES', 'otherTimeRandGMRES','runtimeRandGMRES','sumTimeRandKRPGMRES','otherTimeRandKRPGMRES','runtimeRandKRPGMRES','N','ranksrand','ranksnormal','ranksrandkrp');
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

% if gram
%     exportgraphics(f, 'TTGMRES_Speedup_with_Gram.png')
% else
%     exportgraphics(f, 'TTGMRES_Speedup.png')
%     exportgraphics(f, 'TTGMRES_Speedup.pdf')
% end

f = figure(3);
clf();
f.Position(1:2) = [525,525];
f.Position(3:4) = [525, 355];
ranks = max(ranksrand{1}(1:end-1,:),[],2);

ranks = max(ranksrandkrp{1}(1:end-1,:),[],2);
p3 = plot(ranks, 'r*','markersize',10,'linewidth',2); hold on

p1 = plot(ranks, 'bo','markersize',10,'linewidth',2);
hold on
ranks = max(ranksnormal{1}(1:end-1,:),[],2);
p2 = plot(ranks, 'k+','markersize',10,'linewidth',2);


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

% if gram
%     exportgraphics(f, 'TTGMRES_Ranks_with_Gram.png')
% else
%     exportgraphics(f, 'TTGMRES_Ranks.png')
%     exportgraphics(f, 'TTGMRES_Ranks.pdf')
% end

beep;
pause(1)
beep;


%% Preconditioner
function X = Preconditioner(L,U,p, Y)
    X = Y;
    X{1} = U\(L\Y{1}(p,:));
end
