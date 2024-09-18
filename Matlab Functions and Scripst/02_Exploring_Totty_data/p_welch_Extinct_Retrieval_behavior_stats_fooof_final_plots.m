% Final Plots Granger

% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  01/2024
% Last update: 02/2024


%% Prepare data


%% Load data Extinction

% data organization from mvgc.F_spect --> see data_preprocessing.m
% - Row 1: Baseline freezing
% - Row 2: Baseline not freezing
% - Row 3: CS-TONE freezing
% - Row 4: CS-TONE not freezing
% - Row 5: ITI freezing
% - Row 6: ITI not freezing

for ii = 1:size(pw.over_f.stats_total_power_mean,1)
    pw.over_f.stats_total_power_mean{ii,1}(isempty(pw.over_f.stats_total_power_mean)) = NaN;
    data_2_plot_extinction{ii,1}(:,:,:,ms) = pw.over_f.stats_total_power_mean{ii,1};
end

% Frequency vector
freq_v = pw.over_f.freq_range  ;


%% Load data Retrieval

% data organization from mvgc.F_spect --> see data_preprocessing.m
% - Row 1: Baseline freezing
% - Row 2: Baseline not freezing
% - Row 3: CS-TONE freezing
% - Row 4: CS-TONE not freezing
% - Row 5: ITI freezing
% - Row 6: ITI not freezing
% 
% for ii = 1:size(pw.over_f.stats_total_power_mean,1)
%     pw.over_f.stats_total_power_mean{ii,1}(isempty(pw.over_f.stats_total_power_mean)) = NaN;
%     data_2_plot_retrieval{ii,1}(:,:,:,ms) = pw.over_f.stats_total_power_mean{ii,1};
% end
% 
% % Frequency vector
% freq_v = pw.over_f.freq_range  ;

%%

if ms =
% Mean and SEM
data_2_plot_extinction_mean_SEM = [];

for ii = 1:size(data_2_plot_extinction,1)
    data_2_plot_extinction_mean_SEM{ii,1}  = mean(data_2_plot_extinction{ii,1},4,'omitnan');
    data_2_plot_extinction_mean_SEM{ii,2}  = std(data_2_plot_extinction{ii,1},[],4,'omitnan')./size(data_2_plot_extinction{ii,1},4);
end


% data_2_plot_retrieval_mean_SEM = [];
% 
% for ii = 1:size(data_2_plot_retrieval,1)
%     data_2_plot_retrieval_mean_SEM{ii,1}   = mean(data_2_plot_retrieval{ii,1},4,'omitnan');
%     data_2_plot_retrieval_mean_SEM{ii,2}  = std(data_2_plot_retrieval{ii,1},[],4,'omitnan')./size(data_2_plot_retrieval{ii,1},4);
% end



%% Plot

%data to Plot
data_2_plot = data_2_plot_extinction_mean_SEM;
%data_2_plot = data_2_plot_retrieval_mean_SEM;

figure
set(gcf,'color','w');
% sc = [1,1,960,1200];
% set(gcf, 'Position', sc);
sgtitle({'Welch power spectral density estimate';''});



%Baseline
subplot(1,3,1)
hold on

boundedline(freq_v,data_2_plot{1,1}(1,:),data_2_plot{1,2}(1,:),'LineWidth',2,'color',[.0 .0 .0],'transparency',.4)
boundedline(freq_v,data_2_plot{2,1}(1,:),data_2_plot{2,2}(1,:),'LineWidth',2,'color',[.6 .6 .6],'transparency',.4)
ylabel({'Baseline';'Normalized Power (A.U.)'},'FontWeight','bold','FontSize',12)
xlabel('(Hz)','FontSize',9)
ylim([0 0.1])

subplot(1,3,2)
hold on

boundedline(freq_v,data_2_plot{1,1}(2,:),data_2_plot{1,2}(2,:),'LineWidth',2,'color',[.0 .0 .0],'transparency',.4)
boundedline(freq_v,data_2_plot{2,1}(2,:),data_2_plot{2,2}(2,:),'LineWidth',2,'color',[.6 .6 .6],'transparency',.4)
ylabel({'Baseline';'Normalized Power (A.U.)'},'FontWeight','bold','FontSize',12)
xlabel('(Hz)','FontSize',9)
ylim([0 0.1])

subplot(1,3,3)
hold on

boundedline(freq_v,data_2_plot{1,1}(3,:),data_2_plot{1,2}(3,:),'LineWidth',2,'color',[.0 .0 .0],'transparency',.4)
boundedline(freq_v,data_2_plot{2,1}(3,:),data_2_plot{2,2}(3,:),'LineWidth',2,'color',[.6 .6 .6],'transparency',.4)
ylabel({'Baseline';'Normalized Power (A.U.)'},'FontWeight','bold','FontSize',12)
xlabel('(Hz)','FontSize',9)
ylim([0 0.1])

legend('Freezing','','Non-Freezing','')

%% Save
% newStr1 = id(1:end-8);
% name_1  = strcat(Path,'/',newStr1,'_Granger_Spec_Mean_All_3');
% Path    = files.FilesLoaded{1,1}(ms).folder;
% saveas(gcf,name_1,'png') % save figure
% 
% close all
% 
% clear('ii','jj')
% clear('Combinations_2','Combinations_1','sub_idx','steps','data2plot_mean','cc','session','name_1','newStr1','path','data2plot')
