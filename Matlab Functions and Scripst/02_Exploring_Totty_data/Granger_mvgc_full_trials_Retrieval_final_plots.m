% Final Plots Granger retrieval

% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  01/2024
% Last update: 07/2024


%% Prepare data

%% Load data retrieval

% data organization from mvgc.F_spect --> see data_preprocessing.m
% Row 1: PL  --> IL
% Row 2: IL  --> PL
% Row 3: PL  --> HPC
% Row 4: HPC --> PL
% Row 5: IL  --> HPC
% Row 6: HPC --> IL


data_2_plot_retrieval{1,1}(:,:,ms) = mvgc.stats.Spect_2_12Hz_baseline_2_plot;
data_2_plot_retrieval{2,1}(:,:,ms) = mvgc.stats.Spect_2_12Hz_2_plot_mean;

steps                           = diff(mvgc.parameters.freqs); % according to the fft time window
mvgc.parameters.frex_2_12Hz     = 2:steps(1):12;
mvgc.parameters.frex_idx_2_12Hz = dsearchn(mvgc.parameters.freqs,mvgc.parameters.frex_2_12Hz');

freq_v = mvgc.parameters.freqs(mvgc.parameters.frex_idx_2_12Hz);


%% Mean and SEM


data_2_plot_retrieval_mean_SEM = [];

for ii = 1:size(data_2_plot_retrieval,1)
    data_2_plot_retrieval_mean_SEM{ii,1}  = 10.*mean(data_2_plot_retrieval{ii,1},3,'omitnan');
    data_2_plot_retrieval_mean_SEM{ii,2}  = 10.*std(data_2_plot_retrieval{ii,1},[],3,'omitnan')./size(data_2_plot_retrieval{ii,1},3);
    
end


%% Plot

figure
set(gcf,'color','w');
sc = [1,1,960,1200];
set(gcf, 'Position', sc);
sgtitle({'Pairwise-conditional Granger causality - frequency domain.';''});

subplot(2,2,1)
boundedline(freq_v',data_2_plot_retrieval_mean_SEM{1,1}(3,mvgc.parameters.frex_idx_2_12Hz),data_2_plot_retrieval_mean_SEM{1,2}(3,mvgc.parameters.frex_idx_2_12Hz),'LineWidth',2,'color',[.6 .6 .6],'transparency',.4);
hold on
boundedline(freq_v',data_2_plot_retrieval_mean_SEM{1,1}(4,mvgc.parameters.frex_idx_2_12Hz),data_2_plot_retrieval_mean_SEM{1,2}(4,mvgc.parameters.frex_idx_2_12Hz),'LineWidth',2,'color',[.3 .3 .3],'transparency',.4);
ylabel({'Granger causality'},'FontWeight','bold','FontSize',12)
xlabel('(Hz)','FontSize',12)
title("baseline")
ylim([0 2])

subplot(2,2,2)
boundedline(freq_v',data_2_plot_retrieval_mean_SEM{2,1}(3,mvgc.parameters.frex_idx_2_12Hz),data_2_plot_retrieval_mean_SEM{2,2}(3,mvgc.parameters.frex_idx_2_12Hz),'LineWidth',2,'color',[.6 .6 .6],'transparency',.4);
hold on
boundedline(freq_v',data_2_plot_retrieval_mean_SEM{2,1}(4,mvgc.parameters.frex_idx_2_12Hz),data_2_plot_retrieval_mean_SEM{2,2}(4,mvgc.parameters.frex_idx_2_12Hz),'LineWidth',2,'color',[.3 .3 .3],'transparency',.4);
xlabel('(Hz)','FontSize',12)
title("5 trials")
ylim([0 2])

subplot(2,2,3)
boundedline(freq_v',data_2_plot_retrieval_mean_SEM{1,1}(5,mvgc.parameters.frex_idx_2_12Hz),data_2_plot_retrieval_mean_SEM{1,2}(5,mvgc.parameters.frex_idx_2_12Hz),'LineWidth',2,'color',[.6 .6 .6],'transparency',.4);
hold on
boundedline(freq_v',data_2_plot_retrieval_mean_SEM{1,1}(6,mvgc.parameters.frex_idx_2_12Hz),data_2_plot_retrieval_mean_SEM{1,2}(6,mvgc.parameters.frex_idx_2_12Hz),'LineWidth',2,'color',[.3 .3 .3],'transparency',.4);
ylabel({'Granger causality'},'FontWeight','bold','FontSize',12)
xlabel('(Hz)','FontSize',12)
ylim([0 2])

subplot(2,2,4)
boundedline(freq_v',data_2_plot_retrieval_mean_SEM{2,1}(5,mvgc.parameters.frex_idx_2_12Hz),data_2_plot_retrieval_mean_SEM{2,2}(5,mvgc.parameters.frex_idx_2_12Hz),'LineWidth',2,'color',[.6 .6 .6],'transparency',.4);
hold on
boundedline(freq_v',data_2_plot_retrieval_mean_SEM{2,1}(6,mvgc.parameters.frex_idx_2_12Hz),data_2_plot_retrieval_mean_SEM{2,2}(6,mvgc.parameters.frex_idx_2_12Hz),'LineWidth',2,'color',[.3 .3 .3],'transparency',.4);
xlabel('(Hz)','FontSize',12)
ylim([0 2])

%% Save
% % newStr1 = id(1:end-8);
% % name_1  = strcat(Path,'/',newStr1,'_Granger_Spec_Mean_All_3');
% % Path    = files.FilesLoaded{1,1}(ms).folder;
% % saveas(gcf,name_1,'png') % save figure
% % 
% % close all
% % 
% % clear('ii','jj')
% % clear('Combinations_2','Combinations_1','sub_idx','steps','data2plot_mean','cc','session','name_1','newStr1','path','data2plot')
