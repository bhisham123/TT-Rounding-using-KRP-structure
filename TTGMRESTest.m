%% TTGMRES test
clear;
clc;
maxNumCompThreads(2);
close("all")


%% Add path
addpath("./TTcore/");
addpath("./TTrandomized/");
addpath("./Dataset/")
addpath("./Dataset/nine_cookies/nx_40/");


%% Global parameters

% Number of independent runs
S = 2;

% Number of parameters ranges from 2^2 to 2^l
l = 6; 
N = 2.^(3:l);

% select the cookies data 
dataset = 'new';

%% Load FEM matrices
if strcmp(dataset, 'old')
    list_cookies=[1,2,3,4];
    A{1} = readcoomat('A0.txt');
    A{2} = readcoomat('A1.txt');
    A{3} = readcoomat('A2.txt');
    A{4} = readcoomat('A3.txt');
    A{5} = readcoomat('A4.txt');    
    a0 = readb('b.txt');
    A0 = A{1};
else strcmp(dataset, 'new')
    % selecting cookies out of nine cookies
    list_cookies=[1,2,3,5,7,8,9];
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



%% Run and time experiments
sumTimeNormalGMRES = zeros(length(N),S);
sumTimeRandGMRES = zeros(length(N),S);
sumTimeRandKRPGMRES = zeros(length(N),S);

res_err_normal = zeros(length(N),S);
res_err_rand = zeros(length(N),S);
res_err_randKRP = zeros(length(N),S);


opTimeNormalGMRES = zeros(length(N),S);
opTimeRandGMRES = zeros(length(N),S);
opTimeRandKRPGMRES = zeros(length(N),S);

precTimeNormalGMRES = zeros(length(N),S);
precTimeRandGMRES = zeros(length(N),S);
precTimeRandKRPGMRES = zeros(length(N),S);


remTimeNormalGMRES = zeros(length(N),S);
remTimeRandGMRES = zeros(length(N),S);
remTimeRandKRPGMRES = zeros(length(N),S);


runtimeNormalGMRES = zeros(length(N),S);
runtimeRandGMRES = zeros(length(N),S);
runtimeRandKRPGMRES = zeros(length(N),S);

ranksrand = cell(length(N));
ranksrandkrp = cell(length(N));
ranksnormal = cell(length(N));


for i = 1:length(N)
    n = N(i);

    % Distribution of parameters
    MinD = 1e0;
    MaxD = 1e1;  %for new dataset we use 5

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

    Op             = @(x)       TTsummandsKronOp(C, x);
    deterministic  = @(W,a,tol) TTsum(W, a, tol);
    randomized     = @(W,a,tol) TTsum_Randomize_then_Orthogonalize(W, a, tol);
    randomized_krp = @(W,a,tol) TTsum_Randomize_then_Orthogonalize_KRP(W, a, tol);


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
    end
end

% Discard results of first run
sumTimeNormalGMRES = sumTimeNormalGMRES(:,2:end);
otherTimeNormalGMRES = opTimeNormalGMRES+precTimeNormalGMRES+remTimeNormalGMRES;
otherTimeNormalGMRES = otherTimeNormalGMRES(:,2:end);
runtimeNormalGMRES = runtimeNormalGMRES(:,2:end);



sumTimeRandGMRES = sumTimeRandGMRES(:,2:end);
otherTimeRandGMRES = opTimeRandGMRES+precTimeRandGMRES+remTimeRandGMRES;
otherTimeRandGMRES = otherTimeRandGMRES(:,2:end);
runtimeRandGMRES = runtimeRandGMRES(:,2:end);

sumTimeRandKRPGMRES = sumTimeRandKRPGMRES(:,2:end);
otherTimeRandKRPGMRES = opTimeRandKRPGMRES+precTimeRandKRPGMRES+remTimeRandKRPGMRES;
otherTimeRandKRPGMRES = otherTimeRandKRPGMRES(:,2:end);
runtimeRandKRPGMRES = runtimeRandKRPGMRES(:,2:end);

% file_name = "Cookies_"+dataset+".mat";
% save(file_name, 'sumTimeNormalGMRES', 'otherTimeNormalGMRES', 'runtimeNormalGMRES','sumTimeRandGMRES', 'otherTimeRandGMRES','runtimeRandGMRES','sumTimeRandKRPGMRES','otherTimeRandKRPGMRES','runtimeRandKRPGMRES','N','ranksrand','ranksnormal','ranksrandkrp', 'res_err_normal', 'res_err_rand', 'res_err_randKRP');

%% Plot Results

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
errorbar(N, speedup_krp, neg_krp, pos_krp,'-s','Color','#FF8C00','markersize',10,'linewidth',2); 


%Randomized with TT Structure
speedup = median(sumTimeNormalGMRES, 2) ./ median(sumTimeRandGMRES, 2);
neg = speedup - min(sumTimeNormalGMRES,[],2)./ max(sumTimeRandGMRES,[],2);
pos = max(sumTimeNormalGMRES,[],2)./ min(sumTimeRandGMRES,[],2) - speedup;
errorbar(N, speedup, neg, pos,'b-o','markersize',10,'linewidth',2);
hold on 
speedup = median(runtimeNormalGMRES, 2) ./ median(runtimeRandGMRES, 2);
neg = speedup - min(runtimeNormalGMRES,[],2)./ max(runtimeRandGMRES,[],2);
pos = max(runtimeNormalGMRES,[],2)./ min(runtimeRandGMRES,[],2) - speedup;
errorbar(N, speedup, neg, pos,'-s','Color','#7E2F8E','markersize',10,'linewidth',2);  %'Color','#008000'
set(gca, 'XScale', 'log')
xlabel('Number of parameter samples n_2 = \ldots = n_N (log scale)', 'FontSize', 14)
ylabel('Speedup', 'FontSize', 14)
xticks(N);
set(gca,'FontSize',14)
legend('TT-Sum+Round-KRP','Total speedup-KRP', 'TT-Sum+Round-TTlike','Total speedup-TTlike', 'location', 'best','FontSize', 12.5)
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
legend([p3,p1,p2], 'Randomized TT-GMRES-KRP' , 'Randomized TT-GMRES-TTlike','Naive TT-GMRES', 'location', 'southeast','FontSize', 12.5)
box on
grid on

% file_name = "TTGMRES_"+dataset+".png";
% exportgraphics(f, file_name,'Resolution', 600)



avg_err_normal = mean(res_err_normal,2);
avg_err_rand =  mean(res_err_rand,2);
avg_err_randKRP = mean(res_err_randKRP,2);
f = figure('Position', [100, 100, 520, 400]);
hold on 
plot(N,avg_err_randKRP(1:length(N)),'r*-', 'LineWidth',2,'markersize',10,'DisplayName','Randomized TT-GMRES-KRP'); 
plot(N,avg_err_rand(1:length(N)),'bo-','LineWidth',2,'markersize',10,'DisplayName','Randomized TT-GMRES-TTlike')
plot(N,avg_err_normal(1:length(N)),'k+-','LineWidth',2,'markersize',10,'DisplayName','Naive TT-GMRES')

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

xlabel('Number of parameter samples n_2 = \ldots = n_N (log scale)', 'FontSize', 14)
ylabel('Relative Error', 'FontSize', 14)
xticks(N);
set(gca,'FontSize',14)

legend('location', 'best','FontSize', 12.5)

box on
grid on

% file_name = "TTGMRES_Error_"+dataset+".png";
% exportgraphics(f, file_name,'Resolution', 600)


beep;
pause(1)
beep;

%% Preconditioner
function X = Preconditioner(L,U,p, Y)
    X = Y;
    X{1} = U\(L\Y{1}(p,:));
end
