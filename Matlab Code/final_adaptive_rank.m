%% Synthetic test
clc
clear all ;
close("all")
maxNumCompThreads(1);

%% Add path
addpath("./TTcore/");
addpath("./TTrandomized");
addpath("./Data/");
addpath("./Thermal_Radiation_Transport/")
%% Run and time experiments
S = 3; %number of runs

tols = [1e-1, 1e-2,1e-3,1e-4,1e-5,1e-6,1e-7];

errorsLR = zeros(length(tols),S);
errorsRandLR = zeros(length(tols),S);
errorsRandKRP = zeros(length(tols),S);

timeLR = zeros(length(tols),S);
timeRandLR = zeros(length(tols),S);
timeRandKRP = zeros(length(tols),S);

compLR = zeros(length(tols),S);
compRandLR = zeros(length(tols),S);
compRandKRP =  zeros(length(tols),S);


errorsLR_fix = zeros(length(tols),S);
errorsRandLR_fix = zeros(length(tols),S);
errorsRandOrth_fix = zeros(length(tols),S);
errorsRandKRP_fix = zeros(length(tols),S);

timeLR_fix = zeros(length(tols),S);
timeRandLR_fix = zeros(length(tols),S);
timeRandOrth_fix = zeros(length(tols),S);
timeRandKRP_fix = zeros(length(tols),S);

errorsLR_round = zeros(length(tols),S);
errorsRandLR_round = zeros(length(tols),S);
errorsRandKRP_round = zeros(length(tols),S);

timeLR_round = zeros(length(tols),S);
timeRandLR_round = zeros(length(tols),S);
timeRandKRP_round = zeros(length(tols),S);

compLR_round = zeros(length(tols),S);
compRandLR_round = zeros(length(tols),S);
compRandKRP_round =  zeros(length(tols),S);



seeed = 123;
rng(seeed);

% % Matern data
% d = 100;
file_name = "matern_100_new.mat";
load(file_name)

actual_rank = double(ranks);
I = double(d)*ones(N,1);

TT = cell(N,1);
TT{1} = unfold(tensor(cores.core1),[1,2]);
TT{2} = unfold(tensor(cores.core2),[1,2]);
TT{3} = unfold(tensor(cores.core3),[1,2]);
TT{4} = unfold(tensor(cores.core4),[1,2]);
TT{5} = unfold(tensor(cores.core5),[1,2]);
TT{6} = unfold(tensor(cores.core6),[1,2]);
TT{7} = unfold(tensor(cores.core7),[1,2]);
TT{8} = unfold(tensor(cores.core8),[1,2]);
norm_TT = TTnorm(TT);

fact = 1;
spherical = false;

% filename = ['Matern' '_new.txt'];
% fileID = fopen(filename, 'a');
% 
% r = TTranks(TT);
% fprintf(fileID, 'Ranks:');
% fprintf(fileID, ' %2d', r);
% fprintf(fileID, '\n');

for i=1:length(tols)
    tol = tols(i)
    % fprintf(fileID, 'Tol:%e\n',tol);
    % fprintf(fileID, 'Determinstic TT-Rounding\n');
    for s=1:S
        start = tic;
        Y_Det = TTrounding(TT, tol);
        timeLR(i,s) = toc(start);
        errorsLR(i,s) = TTnorm(TTaxby(1,TT,-1,Y_Det), "OLR") / norm_TT;
        compLR(i,s) =  Compression(TT,Y_Det);

        rankDet = TTranks(Y_Det);
        
        % fprintf(fileID, 'Ranks:');
        % fprintf(fileID, ' %2d', rankDet);
        % fprintf(fileID, '\n');



        % start = tic;
        % Y_Det_fix = TTrounding(TT, 1e-15, rank);
        % timeLR_fix(i,s) = toc(start);
        % errorsLR_fix(i,s) = TTnorm(TTaxby(1,TT,-1,Y_Det), "OLR") / norm_TT;
        % compLR_fix(i,s) =  Compression(TT,Y_Det_fix);
    end
    
    % fprintf(fileID, 'Orthogonalize then Randomize\n');
    for s=1:S
        start = tic;
        Y_RandLR = TTrounding_Orthogonalize_then_Randomize_Adaptive(TT, tol);
        timeRandLR(i,s) = toc(start);
        errorsRandLR(i,s) = TTnorm(TTaxby(1,TT,-1,Y_RandLR), "OLR") / norm_TT;
        compRandLR(i,s) = Compression(TT,Y_RandLR);

        rank = TTranks(Y_RandLR);

        % fprintf(fileID, 'Ranks:');
        % fprintf(fileID, ' %2d', rank);
        % fprintf(fileID, '\n');


        start = tic;
        Y_RandLR_fix = TTrounding_Orthogonalize_then_Randomize(TT, rank);
        timeRandLR_fix(i,s) = toc(start);
        errorsRandLR_fix(i,s) = TTnorm(TTaxby(1,TT,-1,Y_RandLR_fix), "OLR") / norm_TT;

        start = tic;
        Y_RandLR_round = rounding(Y_RandLR,tol);
        timeRandLR_round(i,s) = toc(start);
        errorsRandLR_round(i,s) = TTnorm(TTaxby(1,TT,-1,Y_RandLR_round), "OLR") / norm_TT;
        compRandLR_round(i,s) = Compression(TT,Y_RandLR_round);
    end

    % fprintf(fileID, 'KRP\n');
    for s=1:S
        start = tic;
        Y_RandKtr =  TTroundingKRP_Adaptive(TT, tol);
        timeRandKRP(i,s) = toc(start);
        errorsRandKRP(i,s) = TTnorm(TTaxby(1,TT,-1,Y_RandKtr), "OLR") / norm_TT;
        compRandKRP(i,s) = Compression(TT,Y_RandKtr);

        rank = TTranks(Y_RandKtr);

        % fprintf(fileID, '1 Ranks:');
        % fprintf(fileID, ' %2d', rank);
        % fprintf(fileID, '\n');

        start = tic;
        Y_RandKtr_fix =  TTroundingKRP(TT, rank);
        timeRandKRP_fix(i,s) = toc(start);
        errorsRandKRP_fix(i,s) = TTnorm(TTaxby(1,TT,-1,Y_RandKtr_fix), "OLR") / norm_TT;

        rank = TTranks(Y_RandKtr);
        start = tic;
        Y_RandOrthoLR_fix = TTrounding_Randomize_then_Orthogonalize(TT, rank);
        timeRandOrth_fix(i,s) = toc(start);
        errorsRandOrth_fix(i,s) = TTnorm(TTaxby(1,TT,-1,Y_RandOrthoLR_fix), "OLR") / norm_TT;

        start = tic;
        Y_RandKtr_round =  rounding(Y_RandKtr,tol);
        timeRandKRP_round(i,s) = toc(start);
        errorsRandKRP_round(i,s) = TTnorm(TTaxby(1,TT,-1,Y_RandKtr_round), "OLR") / norm_TT;
        compRandKRP_round(i,s) = Compression(TT,Y_RandKtr_round);
    end    
end


% file_name = "Matern_results_"+num2str(d)+"_"+num2str(S)+"_"+num2str(fact)+"_"+num2str(seeed)+"_"+spherical+"_final.mat";
% 
% save(file_name, 'errorsLR','errorsRandLR','errorsRandKRP','errorsRandLR_round','errorsRandKRP_round','timeLR','timeRandLR','timeRandKRP','timeRandLR_round','timeRandKRP_round','compLR','compRandLR','compRandKRP','compRandLR_round','compRandKRP_round','tols','errorsRandLR_fix','errorsRandOrth_fix','errorsRandKRP_fix','timeRandLR_fix','timeRandOrth_fix','timeRandKRP_fix','S','fact','seeed')


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


%% Vertical subplots

f = figure('Position', [100, 100, 1000, 3500]);


% Compute statistics
med_LR       = median(errorsLR, 2);
min_LR       = min(errorsLR, [], 2);
max_LR       = max(errorsLR, [], 2);

med_RandLR   = median(errorsRandLR, 2);
min_RandLR   = min(errorsRandLR, [], 2);
max_RandLR   = max(errorsRandLR, [], 2);

med_RandKtr  = median(errorsRandKRP, 2);
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
ylim([1e-11, 2e-1])
yticks([1e-11,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1])

xlabel('Tolerance', 'FontSize', 14)
ylabel('Relative Error', 'FontSize', 14)
set(sp1, 'FontSize', 14)
set(sp1, 'YScale', 'log');

lgd1 = legend(sp1, [p2, p3, p1]);
% set(lgd1,'Position',[0.3055 0.8565 0.159 0.067], 'FontSize', 12.5);

box on;
grid on;
% hold off


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
ylabel('Speedup', 'FontSize', 14)
xlabel('Tolerance', 'FontSize', 14)
set(sp2, 'FontSize', 14)
xlim([3,37])

lgd2 = legend(sp2,'show');
% set(lgd2,'Position',[0.3055 0.495929660699357 0.159 0.129356568364611],'FontSize',12.5);
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
% ylim([93 100.5])
xlim([3,37])
lgd3 = legend(sp3,'show');
% set(lgd3,'Position',[0.3055 0.258664245149759 0.159 0.0670241286863271],'FontSize',12.5)

box on;
grid on;
hold off;



%%% plot for Randomized rounding + Compression pass at end
mean_timeRandLR_round = mean(timeRandLR(:,2:end)+timeRandLR_round(:,2:end),2);
mean_timeRandKtr_round = mean(timeRandKRP(:,2:end)+timeRandKRP_round(:,2:end),2);

mean_speedupRandLR_round = mean_timeLR./(mean_timeRandLR_round);
mean_speedupRandKtr_round = mean_timeLR./(mean_timeRandKtr_round);

mean_compRandLR_round = mean(compRandLR_round,2);
mean_compRandKRP_round = mean(compRandKRP_round,2);



% % Create figure and subplot
% f = figure('Position', [100, 100, 1300, 300]);

% Compute statistics
med_LR       = median(errorsLR, 2);
min_LR       = min(errorsLR, [], 2);
max_LR       = max(errorsLR, [], 2);

med_RandLR   = median(errorsRandLR, 2);
min_RandLR   = min(errorsRandLR, [], 2);
max_RandLR   = max(errorsRandLR, [], 2);

med_RandKtr  = median(errorsRandKRP_round, 2);
min_RandKtr  = min(errorsRandKRP_round, [], 2);
max_RandKtr  = max(errorsRandKRP_round, [], 2);

% Compute error bar lengths (top and bottom)
err_low_LR     = med_LR - min_LR;
err_high_LR    = max_LR - med_LR;

err_low_RandLR  = med_RandLR - min_RandLR;
err_high_RandLR = max_RandLR - med_RandLR;

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
ylim([1e-8, 2e-1])
yticks([1e-11,1e-10, 1e-9, 1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1])
xlabel('Tolerance', 'FontSize', 14)
ylabel('Relative Error', 'FontSize', 14)
set(sp4, 'FontSize', 14)
set(sp4, 'YScale', 'log');
lgd4 = legend(sp4,[p2,p3,p1]);
% set(lgd4,'Position',[0.7375 0.859870679466113 0.1675 0.0650134048257373],'FontSize',12);

grid on;
box on;
hold off

sp5 = subplot(3,2,4,'Parent',f);
hold on
plot(5*(1:length(mean_speedupRandKtr)), mean_speedupRandKtr,   'r*-','markersize',10.5,'linewidth',2,'DisplayName', 'TT-KRP (Adap)'); hold on
plot(5*(1:length(mean_speedupRandKtr_round)), mean_speedupRandKtr_round,   's--','Color','#DC143C','markersize',10.5,'linewidth',2,'DisplayName', 'TT-KRP (Adap-R)'); hold on
plot(5*(1:length(mean_speedupRandLR)),  mean_speedupRandLR,'b+-','markersize',10.5,'linewidth',2,'DisplayName', 'TT-Orth-Rand (Adap)')
% plot(5*(1:length(mean_speedupRandLR_round)),  mean_speedupRandLR_round,'x--','Color','#0047AB','markersize',10.5,'linewidth',2,'DisplayName', 'TT-Orth-Rand (Adap-R)')
plot(5*(1:length(mean_speedupLR)), mean_speedupLR,'ko-','markersize',10.5,'linewidth',2,'DisplayName', 'TT-Rounding')

xticks(5*(1:7))
xticklabels({'1e-1','1e-2', '1e-3', '1e-4', '1e-5','1e-6','1e-7'})
ylabel('Speedup', 'FontSize', 14)
xlabel('Tolerance', 'FontSize', 14)
set(sp5, 'FontSize', 14)
xlim([3,37])
lgd5=legend(sp5, 'show');
% set(lgd5,'Position',[0.7325 0.516707140592119 0.1725 0.10857908847185],'FontSize',12.5);
grid on
box on;
hold off

sp6 = subplot(3,2,6,'Parent',f);
p1 = plot(5*(1:length(mean_compLR)), mean_compLR,'ko-','markersize',10.5,'linewidth',2,'DisplayName', 'TT-Rounding');hold on
p2 = plot(5*(1:length(mean_compRandKRP_round)), mean_compRandKRP_round,   's','Color','#DC143C','markersize',10.5,'linewidth',2,'DisplayName', 'TT-KRP (Adap-R)'); 
p3 = plot(5*(1:length(mean_compRandLR_round)),  mean_compRandLR_round,'x','Color','#0047AB','markersize',10.5,'linewidth',2,'DisplayName', 'TT-Orth-Rand (Adap-R)');
% p3 = plot(5*(1:length(mean_compRandLR)),  mean_compRandLR,'b+-','markersize',10.5,'linewidth',2,'DisplayName', 'TT-Orth-Rand (Adap)');

xticks(5*(1:7))
xticklabels({'1e-1','1e-2', '1e-3', '1e-4', '1e-5','1e-6','1e-7'})
yticks([1e0, 1e1,1e2,1e3,1e4,1e5])

ylabel('Compression', 'FontSize', 14)
xlabel('Tolerance', 'FontSize', 14)
set(sp6, 'FontSize', 14)
set(sp6, 'YScale', 'log')

xlim([3,37])
lgd6 = legend(sp6,[p2,p3,p1]);
% set(lgd6,'Position',[0.7375 0.260674969010349 0.1675 0.0650134048257373],'FontSize',12);

grid on
box on 
hold off
    
function ratio = Compression(X,Y)
    [~,d,r] = TTsizes(X);
    [~,d,l] = TTsizes(Y);
    % Compute original and rounded sizes
    original_size = sum(r(1:end-1) .* d .* r(2:end));
    rounded_size = sum(l(1:end-1) .* d .* l(2:end));
    ratio = original_size/rounded_size;
end

function X = rounding(X,tol)
    % Cores of X are left orthonromal except the last one
    [N,I,~] = TTsizes(X);
    normX = norm(X{N}, 'fro');
    tau = tol * normX / sqrt(N - 1);                    % Truncation threshold
    for n = N:-1:2
        [Q,R] = qr(v2h(X{n}, I(n))', 0);
        [U, S, V] = svd(R, 'econ');
        rank = trunc(diag(S), tau);
        X{n} = h2v( (Q*U(:, 1:rank))', I(n));
        try
        X{n-1} = X{n-1} * (V(:, 1:rank) * S(1:rank, 1:rank));
        catch
            warning('Some problem');
            break;
        end
    end
end


   
    
    