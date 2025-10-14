clear all ;
close("all")
maxNumCompThreads(1);

%% Add path
addpath("./TTcore/");
addpath("./TTrandomized");

%% accuracy experiments

Runs = 1000; %number of runs
Samples = 10:10:200;

Colors = {[0 0.4470 0.7410],  [0.4660 0.6740 0.1880],[0.8500 0.3250 0.0980], [0.4940 0.1840 0.5560]};

legend_handles = gobjects(4,1);

f = figure();
hold on
legend_handles(1) = yline(1, 'k--', 'LineWidth', 2, 'DisplayName', 'Actual Norm');   % Reference line at y=1


est =cell(3,1);

k = 1;
for N = 15:-5:5
    N
    errors = zeros(length(Samples),Runs);
    NormEst = zeros(length(Samples),Runs);

    d = 100;  % Size of each mode
 
    perturbation = 1e-5;
    
    % Rank
    r1 = 50;  
    
    % Perturbation Rank
    r2 = 50; 

    % % Mode dimensions and ranks
    I = d*ones(N,1);
    R1 = [1; r1*ones(N-1,1); 1];
    R2 = [1; r2*ones(N-1,1); 1];
    TT1 = TTrand(I,R1);
    TT2 = TTrand(I,R2);
    
    norm_TT1 = TTnorm(TT1);
    norm_TT2 = TTnorm(TT2);
    
    TT = TTaxby(1/norm_TT1,TT1,perturbation/norm_TT2,TT2);
    norm_TT = TTnorm(TT);
    
    
    for i = 1: length(Samples)
        l = Samples(i);
        for j = 1:Runs
            W = KRPPartialContractionsRL(TT,l);
            S = TT{1}*W{1};
            NormEst(i,j) = norm(S,'fro')/sqrt(l);
            errors(i,j) = abs(norm_TT-NormEst(i,j))/norm_TT;
        end
    end

    
    est{k}     = mean(NormEst, 2);
    % Prepare shaded region
    lower = min(NormEst, [], 2);
    upper = max(NormEst, [], 2);

    x = [Samples'; flipud(Samples')];
    y = [lower; flipud(upper)];

    % Plot shaded region without legend entry
    fill(x, y, Colors{k}, 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    k = k+1
end

legend_handles(2) = plot(Samples, est{1}, 'Color', Colors{1}, 'LineWidth', 2, 'DisplayName', "d = " + num2str(15));
legend_handles(3) = plot(Samples, est{2}, 'Color', Colors{2}, 'LineWidth', 2, 'DisplayName', "d = " + num2str(10));
legend_handles(4) = plot(Samples, est{3}, 'Color', Colors{3}, 'LineWidth', 2, 'DisplayName', "d = " + num2str(5));

set(gca, 'YScale', 'log')
ylabel('Estimated Norm')
xlabel('Number of Samples')
set(gca, 'FontSize', 16)
xlim([10,200])

legend_handles = flipud(legend_handles);
labels = arrayfun(@(h) get(h, 'DisplayName'), legend_handles, 'UniformOutput', false);
legend(legend_handles, labels,'Location', 'southEast');
box on

% file_name = "Norm_estimte.png";
% exportgraphics(f, file_name,'Resolution', 600);



