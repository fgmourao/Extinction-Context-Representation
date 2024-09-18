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

for ii = 1:size(pw.behavior.stats.spectrum_2_12Hz_mean,1)
    pw.stats_total_power_mean{ii,2}(isempty(pw.behavior.stats.spectrum_2_12Hz_mean)) = NaN;
    data_2_plot_extinction{ii,1}(:,:,:,ms) = pw.behavior.stats.spectrum_2_12Hz_mean{ii,2};
end

% Frequency vector
freq_v = pw.behavior.parameters.freq_(pw.behavior.parameters.frex_idx_2_12Hz);


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
%     pw.over_f.stats_total_power_mean{ii,2}(isempty(pw.over_f.stats_total_power_mean)) = NaN;
%     data_2_plot_retrieval{ii,1}(:,:,:,ms) = pw.over_f.stats_total_power_mean{ii,2};
% end
% 
% % Frequency vector
% freq_v = pw.behavior.parameters.freq_(pw.behavior.parameters.frex_idx_2_12Hz);

%%

if ms == length(files.FilesLoaded{1, 1})

% Mean and SEM

data_2_plot_extinction_mean_SEM = [];

for ii = 1:size(data_2_plot_extinction,1)
    data_2_plot_extinction_mean_SEM{ii,1}  = 10.*mean(data_2_plot_extinction{ii,1},4,'omitnan');
    data_2_plot_extinction_mean_SEM{ii,2}  = 10.*std(data_2_plot_extinction{ii,1},[],4,'omitnan')./size(data_2_plot_extinction{ii,1},4);
    
    %SEM shades
    curve1_ext{ii,1} = data_2_plot_extinction_mean_SEM{ii, 1} + data_2_plot_extinction_mean_SEM{ii, 2};
    curve2_ext{ii,1} = data_2_plot_extinction_mean_SEM{ii, 1} - data_2_plot_extinction_mean_SEM{ii, 2};
    data_2_plot_extinction_mean_SEM{ii,3} = cat(2,curve1_ext{ii,1},flip(curve2_ext{ii,1},2));
    data_2_plot_extinction_mean_SEM{ii,4} = [freq_v', fliplr(freq_v')];

end


% % data_2_plot_retrieval_mean_SEM = [];
% % 
% % for ii = 1:size(data_2_plot_retrieval,1)
% %     data_2_plot_retrieval_mean_SEM{ii,1}   = 10.*mean(data_2_plot_retrieval{ii,1},4,'omitnan');
% %     data_2_plot_retrieval_mean_SEM{ii,2}  = 10.*std(data_2_plot_retrieval{ii,1},[],4,'omitnan')./size(data_2_plot_retrieval{ii,1},4);
% % 
% %     
% %     % SEM shades
% %     curve1_ext{ii,1} = data_2_plot_retrieval_mean_SEM{ii, 1} + data_2_plot_retrieval_mean_SEM{ii, 2};
% %     curve2_ext{ii,1} = data_2_plot_retrieval_mean_SEM{ii, 1} - data_2_plot_retrieval_mean_SEM{ii, 2};
% %     data_2_plot_retrieval_mean_SEM{ii,3} = cat(2,curve1_ext{ii,1},flip(curve2_ext{ii,1},2));
% %     data_2_plot_retrieval_mean_SEM{ii,4} = [freq_v', fliplr(freq_v')];
% % 
% % end
% 
% 

%% Plot

%data to Plot
data_2_plot = data_2_plot_extinction_mean_SEM;
%data_2_plot = data_2_plot_retrieval_mean_SEM;

figure
set(gcf,'color','w');
sc = [1,1,960,1200];
set(gcf, 'Position', sc);
sgtitle({'Welch power spectral density estimate';''});



%Baseline
subplot(3,3,1)
hold on

fill(data_2_plot{1,4},data_2_plot{1,3}(1,:),'k','FaceAlpha',0.5,'EdgeColor','none')
plot(freq_v,data_2_plot{1,1}(1,:),'color','k','linew',2)
fill(data_2_plot{2,4},data_2_plot{2,3}(1,:),[.4 .4 .4],'FaceAlpha',0.3,'EdgeColor','none')
plot(freq_v,data_2_plot{2,1}(1,:),'color','k','linew',2)
ylabel({'Baseline';'Normalized Power (A.U.)'},'FontWeight','bold','FontSize',12)
xlabel('(Hz)','FontSize',9)
ylim([0 0.12])

subplot(3,3,2)
hold on

fill(data_2_plot{1,4},data_2_plot{1,3}(2,:),'k','FaceAlpha',0.5,'EdgeColor','none')
plot(freq_v,data_2_plot{1,1}(2,:),'color','k','linew',2)
fill(data_2_plot{2,4},data_2_plot{2,3}(2,:),[.4 .4 .4],'FaceAlpha',0.3,'EdgeColor','none')
plot(freq_v,data_2_plot{2,1}(2,:),'color',[.4 .4 .4],'linew',2)
xlabel('(Hz)','FontSize',9)
ylim([0 0.12])

subplot(3,3,3)
hold on

fill(data_2_plot{1,4},data_2_plot{1,3}(3,:),'k','FaceAlpha',0.5,'EdgeColor','none')
plot(freq_v,data_2_plot{1,1}(3,:),'color','k','linew',2)
fill(data_2_plot{2,4},data_2_plot{2,3}(3,:),[.4 .4 .4],'FaceAlpha',0.3,'EdgeColor','none')
plot(freq_v,data_2_plot{2,1}(3,:),'color',[.4 .4 .4],'linew',.2)
xlabel('(Hz)','FontSize',9)
ylim([0 0.12])

legend('Freezing','','Non-Freezing','')



%CS-TONE
subplot(3,3,1+3)
hold on

fill(data_2_plot{3,4},data_2_plot{3,3}(1,:),[0.6350, 0.0780, 0.1840],'FaceAlpha',0.2,'EdgeColor','none')
plot(freq_v,data_2_plot{3,1}(1,:),'color',[0.6350, 0.0780, 0.1840],'linew',2)
fill(data_2_plot{4,4},data_2_plot{4,3}(1,:),[1, .4, .2],'FaceAlpha',0.3,'EdgeColor','none')
plot(freq_v,data_2_plot{4,1}(1,:),'color',[1, .4, .2],'linew',2)
ylabel({'CS-TONE';'Normalized Power (A.U.)'},'FontWeight','bold','FontSize',12)
xlabel('(Hz)','FontSize',9)
ylim([0 0.12])

subplot(3,3,2+3)
hold on

fill(data_2_plot{3,4},data_2_plot{3,3}(2,:),[0.6350, 0.0780, 0.1840],'FaceAlpha',0.2,'EdgeColor','none')
plot(freq_v,data_2_plot{3,1}(2,:),'color',[0.6350, 0.0780, 0.1840],'linew',2)
fill(data_2_plot{4,4},data_2_plot{4,3}(2,:),[1, .4, .2],'FaceAlpha',0.3,'EdgeColor','none')
plot(freq_v,data_2_plot{4,1}(2,:),'color',[1, .4, .2],'linew',2)
xlabel('(Hz)','FontSize',9)
ylim([0 0.12])

subplot(3,3,3+3)
hold on

fill(data_2_plot{3,4},data_2_plot{3,3}(3,:),[0.6350, 0.0780, 0.1840],'FaceAlpha',0.2,'EdgeColor','none')
plot(freq_v,data_2_plot{3,1}(3,:),'color',[0.6350, 0.0780, 0.1840],'linew',2)
fill(data_2_plot{4,4},data_2_plot{4,3}(3,:),[1, .4, .2],'FaceAlpha',0.3,'EdgeColor','none')
plot(freq_v,data_2_plot{4,1}(3,:),'color',[1, .4, .2],'linew',2)
xlabel('(Hz)','FontSize',9)
ylim([0 0.12])
legend('Freezing','','Non-Freezing','')

%ITI
subplot(3,3,1+6)
hold on

fill(data_2_plot{5,4},data_2_plot{5,3}(1,:),[.2, .2, 1],'FaceAlpha',0.5,'EdgeColor','none')
plot(freq_v,data_2_plot{5,1}(1,:),'color',[.2, .2, 1],'linew',2)
fill(data_2_plot{6,4},data_2_plot{6,3}(1,:),[.2, .2, 1],'FaceAlpha',0.2,'EdgeColor','none')
plot(freq_v,data_2_plot{6,1}(1,:),'color',[.2, .2, 1],'linew',2)
ylabel({'ITI';'Normalized Power (A.U.)'},'FontWeight','bold','FontSize',12)
xlabel('(Hz)','FontSize',9)
ylim([0 0.12])

subplot(3,3,2+6)
hold on

fill(data_2_plot{5,4},data_2_plot{5,3}(2,:),[.2, .2, 1],'FaceAlpha',0.5,'EdgeColor','none')
plot(freq_v,data_2_plot{5,1}(2,:),'color',[.2, .2, 1],'linew',2)
fill(data_2_plot{6,4},data_2_plot{6,3}(2,:),[.2, .2, 1],'FaceAlpha',0.2,'EdgeColor','none')
plot(freq_v,data_2_plot{6,1}(2,:),'color',[.2, .2, 1],'linew',2)
xlabel('(Hz)','FontSize',9)
ylim([0 0.12])

subplot(3,3,3+6)
hold on

fill(data_2_plot{5,4},data_2_plot{5,3}(3,:),[.2, .2, 1],'FaceAlpha',0.5,'EdgeColor','none')
plot(freq_v,data_2_plot{5,1}(3,:),'color',[.2, .2, 1],'linew',2)
fill(data_2_plot{6,4},data_2_plot{6,3}(3,:),[.2, .2, 1],'FaceAlpha',0.2,'EdgeColor','none')
plot(freq_v,data_2_plot{6,1}(3,:),'color',[.2, .2, 1],'linew',2)
xlabel('(Hz)','FontSize',9)
ylim([0 0.12])

legend('Freezing','','Non-Freezing','')

end

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
