%% Synthetic test
clear all ;
close("all")
maxNumCompThreads(1);

%% Add path
addpath("./TTcore/");
addpath("./TTrandomized");
addpath("./Data/");





%% Run and time experiments
S = 1; %number of runs

tol = 1e-1;


errorsLR = zeros(1,S);
errorsRandLR = zeros(1,S);
errorsRandKtr = zeros(1,S);

timeLR = zeros(1,S);
timeRandLR = zeros(1,S);
timeRandKtr = zeros(1,S);




% % % Matern data
% d = 100;
% file_name = "matern_data.mat";
% load(file_name)
% 
% actual_rank = double(ranks);
% I = double(d)*ones(N,1);
% 
% TT = cell(N,1);
% TT{1} = unfold(tensor(cores.core1),[1,2]);
% TT{2} = unfold(tensor(cores.core2),[1,2]);
% TT{3} = unfold(tensor(cores.core3),[1,2]);
% TT{4} = unfold(tensor(cores.core4),[1,2]);
% TT{5} = unfold(tensor(cores.core5),[1,2]);
% TT{6} = unfold(tensor(cores.core6),[1,2]);
% TT{7} = unfold(tensor(cores.core7),[1,2]);
% TT{8} = unfold(tensor(cores.core8),[1,2]);
% norm_TT = TTnorm(TT);


%%% H2O
d = 16;
file_name = "h2o_stog3_tol1e16.mat";
load(file_name);


rank = double(ranks);
I = double(d)*ones(N,1);
TT = cell(N,1);
TT{1} = unfold(tensor(cores{1}),[1,2]);
TT{2} = unfold(tensor(cores{2}),[1,2]);
TT{3} = unfold(tensor(cores{3}),[1,2]);
TT{4} = unfold(tensor(cores{4}),[1,2]);
TT{5} = unfold(tensor(cores{5}),[1,2]);
TT{6} = unfold(tensor(cores{6}),[1,2]);
TT{7} = unfold(tensor(cores{7}),[1,2]);
norm_TT = TTnorm(TT);


% % % %%% N_2 molecule
% d = 16;
% file_name = "n2_tol1e15.mat";
% load(file_name);
% 
% rank = double(ranks);
% I = double(d)*ones(N,1);
% TT = cell(N,1);
% TT{1} = unfold(tensor(cores{1}),[1,2]);
% TT{2} = unfold(tensor(cores{2}),[1,2]);
% TT{3} = unfold(tensor(cores{3}),[1,2]);
% TT{4} = unfold(tensor(cores{4}),[1,2]);
% TT{5} = unfold(tensor(cores{5}),[1,2]);
% TT{6} = unfold(tensor(cores{6}),[1,2]);
% TT{7} = unfold(tensor(cores{7}),[1,2]);
% TT{8} = unfold(tensor(cores{8}),[1,2]);
% TT{9} = unfold(tensor(cores{9}),[1,2]);
% TT{10} = unfold(tensor(cores{10}),[1,2]); 
% norm_TT = TTnorm(TT);



for i=1:1
    for s=1:S
        start = tic;
        Y_Det = TTrounding(TT, tol);
        timeLR(i,s) = toc(start);
        errorsLR(i,s) = TTnorm(TTaxby(1,TT,-1,Y_Det), "OLR") / norm_TT;
    end
    
    for s=1:S
        start = tic;
        Y_RandLR = TTrounding_Orthogonalize_then_Randomize_Adaptive(TT, tol);
        timeRandLR(i,s) = toc(start);
        errorsRandLR(i,s) = TTnorm(TTaxby(1,TT,-1,Y_RandLR), "OLR") / norm_TT;
    end
    
    for s=1:S
        start = tic;
        % Y_RandKtr =  TTrounding_KTR_Adaptive_Exact(TT, tol);
        Y_RandKtr =  TTrounding_KTR_Adaptive_Est(TT,norm_TT, tol);
        timeRandKtr(i,s) = toc(start);
        errorsRandKtr(i,s) = TTnorm(TTaxby(1,TT,-1,Y_RandKtr), "OLR") / norm_TT;
    end
        
end

    

mean_errLR = mean(errorsLR,2);
mean_errRandLR = mean(errorsRandLR,2);
mean_errRandKtr = mean(errorsRandKtr,2);



 

f = figure();
subplot(1,2,1,'Position', [0.1, 0.3, 0.35, 0.6]);
    
% Data
Err = [mean_errLR, mean_errRandLR, mean_errRandKtr];
    
% Plot the bars separately for each color
h1 = bar(1, Err(1));  % Plot the first bar at position 1
hold on;
h2 = bar(2, Err(2));  % Plot the second bar at position 2
h3 = bar(3, Err(3));  % Plot the second bar at position 2

% Set colors for each bar
h1.FaceColor = [0, 0, 0];  % Black for the first bar
h2.FaceColor = [0, 0, 1];  % Blue for the first bar
h3.FaceColor = [1, 0, 0];  % Red for the second bar

% Set the Y-axis to logarithmic scale
ax = gca;
ax.YScale = 'log';
ylim([min(Err)/2, max(Err)*2]);

% Set the labels and legend
ylabel('Mean Relative Error', 'FontSize', 18);
legend([h1, h2, h3], {'TT-Rounding','TT-Orth-Rand', 'Rand-KRP'}, 'FontSize', 8);

% Remove extra ticks
set(gca, 'XTick', 1:2);
set(gca, 'XTickLabel', {'', ''});  % Optionally remove the X-axis labels


legend('location','northwest','FontSize',15);
set(gca,'FontSize',18)
  

    
 
mean_timeLR = mean(timeLR,2);
mean_speedupLR = mean_timeLR./mean(timeLR,2);
mean_speedupRandLR = mean(timeLR,2)./mean(timeRandLR,2);
mean_speedupRandKtr = mean_timeLR./mean(timeRandKtr,2);
   


    
subplot(2,2,2,'Position', [0.55, 0.3, 0.35, 0.6]);

plot(1:length(mean_speedupLR), mean_speedupLR','ko-','markersize',12,'linewidth',2)
hold on
plot(1:length(mean_speedupRandLR), mean_speedupRandLR','bs-','markersize',12,'linewidth',2)
plot(1:length(mean_speedupRandKtr), mean_speedupRandKtr',   'r*-' ,'markersize',12,'linewidth',2)
hold off
ylabel('Speedup', 'FontSize', 18)
xticks(1:length(mean_errLR));  % Specify the positions of the ticks
set(gca,'FontSize',18)
sgtitle(['Tol = ',num2str(tol)],'fontsize',18);
% Remove x-axis ticks
xticks([]);
    
   
    
    