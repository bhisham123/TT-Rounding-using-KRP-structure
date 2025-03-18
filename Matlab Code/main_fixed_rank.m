%% Synthetic test
clear all ;
close("all")
maxNumCompThreads(1);

%% Add path
addpath("../TTcore/");
addpath("../TTrandomized");
addpath("../Data/");




%% Run and time experiments
S = 2; %number of runs

tol = 1e-2;
over_sample = [3,5,7]; %[1,2,3,5,7];
L = length(over_sample);

errorsLR = zeros(L,S);
errorsRandLR = zeros(L,S);
errorsRandOrth = zeros(L,S);
errorsTwoSide = zeros(L,S);
errorsRandKtr = zeros(L,S);
errorsTwoSideKtr = zeros(L,S);

timeLR = zeros(L,S);
timeRandLR = zeros(L,S);
timeRandOrth = zeros(L,S);
timeTwoSide = zeros(L,S);
timeRandKtr = zeros(L,S);
timeTwoSideKtr = zeros(L,S);



% % % Matern data
d = 100;
file_name = "matern_data.mat";
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


%%% H2O
% d = 16;
% file_name = "h2o_stog3_tol1e16.mat";
% load(file_name);
% 
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
% norm_TT = TTnorm(TT);
% 
% 
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
% 
% 
 
X = TTrounding(TT, tol);
for i=1:L

    p = over_sample(i)
    
    rank(1) = 1;
    rankp(1) = 1;
    for  j = 1:length(X)-1
            rank(j+1) = size(X{j},2);
         rankp(j+1) = size(X{j},2)+ p;
    end
    rank(length(X)+1) = 1;
    rankp(length(X)+1) = 1;

    
    for s=1:S
        start = tic;
        Y_Det = TTrounding(TT, 1e-15, rank);
        timeLR(i,s) = toc(start);
        errorsLR(i,s) = TTnorm(TTaxby(1,TT,-1,Y_Det), "OLR") / norm_TT;
    end
    
    for s=1:S
        start = tic;
        Y_OrthRand = TTrounding_Orthogonalize_then_Randomize(TT, rankp);
        timeRandLR(i,s) = toc(start);
        errorsRandLR(i,s) = TTnorm(TTaxby(1,TT,-1,Y_OrthRand), "OLR") / norm_TT;
    end
    
    for s=1:S
        start = tic;
        Y_RandOrth = TTrounding_Randomize_then_Orthogonalize(TT, rankp);
        timeRandOrth(i,s) = toc(start);
        errorsRandOrth(i,s) = TTnorm(TTaxby(1,TT,-1,Y_RandOrth), "OLR") / norm_TT;
    end
    
    for s=1:S
        start = tic;
        Y_TwoSide = TTrounding_Two_Sided_Randomization(TT, rankp);
        timeTwoSide(i,s) = toc(start);
        errorsTwoSide(i,s) = TTnorm(TTaxby(1,TT,-1,Y_TwoSide), "OLR") / norm_TT;
    end
    
    for s=1:S
        start = tic;
        Y_RandKtr =  TTrounding_KTR(TT, rankp);
        timeRandKtr(i,s) = toc(start);
        errorsRandKtr(i,s) = TTnorm(TTaxby(1,TT,-1,Y_RandKtr), "OLR") / norm_TT;
    end


    for s=1:S
        start = tic;
        Y_TwoSide = TTrounding_Two_Sided_Randomization_KRP(TT, rankp);
        timeTwoSideKtr(i,s) = toc(start);
        errorsTwoSideKtr(i,s) = TTnorm(TTaxby(1,TT,-1,Y_TwoSide), "OLR") / norm_TT;
    end
        
end



mean_errLR = mean(errorsLR,2);
mean_errRandLR = mean(errorsRandLR,2);
mean_errRandOrth = mean(errorsRandOrth,2);
mean_errTwoSide = mean(errorsTwoSide,2);
mean_errRandKtr = mean(errorsRandKtr,2);
mean_errTwoSideKtr = mean(errorsTwoSideKtr,2);

    
f = figure();
subplot(1,2,1,'Position', [0.1, 0.3, 0.35, 0.6]);
    
Err = [mean_errLR,mean_errRandLR,mean_errRandOrth,mean_errRandKtr,mean_errTwoSide,mean_errTwoSideKtr];
h = bar(1:length(mean_errLR),Err);
ax = gca;
ax.YScale = 'log';

colors = {'black','blue','#993399','#FF3300','#99CC33', "#4DBEEE"};
% Set the colors for each bar group
for k = 1:length(h)
    h(k).FaceColor = colors{k};  % Set the color for each bar group
end

% title(['tol = ',num2str(tol)])
xlabel('Over Sampling (p)', 'FontSize', 18)
ylabel('Mean Relative Error', 'FontSize', 18)
legend('TT-Rounding','Orth-then-Rand', 'Rand-TTlike','Rand-KRP','Two-Sided-TTlike', 'Two-Sided-KRP','FontSize', 16)
legend('location','northeast');
set(gca,'FontSize',18)
xticks(1:length(mean_errLR));  % Specify the positions of the ticks
xticklabels(over_sample);  % Set the labels
ylim([min(Err(:))*0.85,max(Err(:))*1.2])
        

    

mean_timeLR = mean(timeLR,2);
mean_speedupLR = mean_timeLR./mean(timeLR,2);
mean_speedupRandLR = mean_timeLR./mean(timeRandLR,2);
mean_speedupRandOrth = mean_timeLR./mean(timeRandOrth,2);
mean_speedupTwoSide = mean_timeLR./mean(timeTwoSide,2);
mean_speedupRandKtr = mean_timeLR./mean(timeRandKtr,2);
mean_speedupTwoSideKtr = mean_timeLR./mean(timeTwoSideKtr,2);


  
subplot(2,2,2,'Position', [0.55, 0.3, 0.35, 0.6]);

plot(1:length(mean_speedupLR), mean_speedupLR', 'o-', 'Color', 'black', 'MarkerSize', 12, 'LineWidth', 3)
hold on

plot(1:length(mean_speedupRandLR), mean_speedupRandLR', 's-','color', 'blue','markersize',12,'linewidth',3)
plot(1:length(mean_speedupRandOrth), mean_speedupRandOrth',  'd-','color','#993399' ,'markersize',12,'linewidth',3)
plot(1:length(mean_speedupRandKtr), mean_speedupRandKtr',   '*-' ,'color', '#FF3300','markersize',12,'linewidth',3)
plot(1:length(mean_speedupTwoSide),mean_speedupTwoSide',  '+-', 'color','#99CC33' ,'markersize',12,'linewidth',3)
plot(1:length(mean_speedupTwoSideKtr),mean_speedupTwoSideKtr',  'x-' ,'color', "#4DBEEE",'markersize',12,'linewidth',3)
hold off
xlabel('Over Sampling (p)', 'FontSize', 18)
ylabel('Speedup', 'FontSize', 18)

set(gca,'FontSize',18)
xlim([0.8, length(over_sample)+0.2])
xticks(1:length(mean_errLR));  % Specify the positions of the ticks
xticklabels(over_sample); 


sgtitle(['Tol = ',num2str(tol)],'FontSize', 18);  %', b = ',num2str(b)




