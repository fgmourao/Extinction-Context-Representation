% Final Plots Granger

% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  01/2024
% Last update: 02/2024


%% Prepare data


%% Loading data Extinction
%data_2_plot_extinction(:,:,:,ms) = mvgc.stats.Spect_2_12Hz_2_plot_mean;

%% Load data Retrieval

% data organization from mvgc.F_spect --> see data_preprocessing.m
% - Row 1: Baseline freezing
% - Row 2: Baseline not freezing
% - Row 3: CS-TONE freezing
% - Row 4: CS-TONE not freezing
% - Row 5: ITI freezing
% - Row 6: ITI not freezing

for ii = 1:size(mvgc.stats.spectrum_2_12Hz_mean,1)
    mvgc.stats.spectrum_2_12Hz_mean{ii,1}(isempty(mvgc.stats.spectrum_2_12Hz_mean{ii,1})) = NaN;
    data_2_plot{ii,1}(:,:,ms) = mvgc.stats.spectrum_2_12Hz_mean{ii,1};
end


%% Frequency vector
% 
% freq_v = mvgc.parameters.freqs(mvgc.parameters.frex_idx_2_12Hz);
% 
% %% Mean and SEM
% 
% for ii = 1:size(data_2_plot{1,1},1)
%     data_2_plot_mean{ii,1} = mean(data_2_plot{ii,1},3,'omitnan');
%     data_2_plot_SEM{ii,1}  = std(data_2_plot{ii,1},[],3,'omitnan')./size(data_2_plot{ii,1},3);
% 
%     % SEM shades
%     curve1{ii,1} = squeeze(plot_1_mean(1,:) + plot_1_SEM(1,:));
%     curve2{ii,1} = squeeze(plot_1_mean(1,:) - plot_1_SEM(1,:));
%     x2{ii,1} = [freq_v', fliplr(freq_v')];
%     inBetween{ii,1} = [curve1{ii,1}, fliplr(curve2{ii,1})];
% 
% end
% 
% %% Choose data toplot
% 
% % idx = 1 Figure --> within CS
% % idx = 1 Figure --> outside CS
% 
% for ii = 1:2
%     % Figure 1 - within CS
%     % Figure 2 - outside CS
% 
%     plot_1_mean = data_2_plot_mean{1,1};
%     plot_1_SEM  = data_2_plot_SEM{1,1};
% 
%     plot_2_mean = data_2_plot_mean{2,ii};
%     plot_2_SEM  = data_2_plot_SEM{2,ii};
% 
%     plot_3_mean = data_2_plot_mean{3,ii};
%     plot_3_SEM  = data_2_plot_SEM{3,ii};
% 
% 
%     figure(ii+4)
%     set(gcf,'color','w');
%     set(gcf, 'Position', get(0, 'Screensize'));
% 
%     if ii == 1
%         sgtitle({'Multivariate Granger Causality CS-Tone - VAR model estimation regression mode (LWR). Order ~= 35';[]})
%     else
%         sgtitle({'Multivariate Granger Causality Outside CS-Tone - VAR model estimation regression mode (LWR). Order ~= 35';[]})
% 
%     end
% 
%     % Baseline
%     subplot(3,3,1)
% 
% 
%     fill(x2, inBetween, 'k','FaceAlpha',0.6,'EdgeColor','none');
%     hold on
%     plot(mvgc.parameters.freqs(1:mvgc.parameters.fres),squeeze(plot_1_mean(1,:)),'linew', 2,'Color','k')
% 
%     curve1 = squeeze(plot_1_mean(2,:) + plot_1_SEM(2,:));
%     curve2 = squeeze(plot_1_mean(2,:) - plot_1_SEM(2,:));
%     x2 = [mvgc.parameters.freqs(1:mvgc.parameters.fres)', fliplr(mvgc.parameters.freqs(1:mvgc.parameters.fres)')];
%     inBetween = [curve1, fliplr(curve2)];
%     fill(x2, inBetween, 'k','FaceAlpha',0.2,'EdgeColor','none');
%     plot(mvgc.parameters.freqs(1:mvgc.parameters.fres),squeeze(plot_1_mean(2,:)),'linew', 2,'Color',[ 0 0 0 .5])
%     xline(6)
%     xlim([2 15])
%     ylim([0 .40])
%     xlabel('(Hz)','FontSize',9), ylabel({'Baseline'; 'Granger Prediction'},'FontSize',11)
%     legend('','PL --> IL','','IL --> PL')
%     title('PL <--> IL')
%     box off
% 
%     subplot(3,3,2)
%     curve1 = squeeze(plot_1_mean(3,:) + plot_1_SEM(3,:));
%     curve2 = squeeze(plot_1_mean(3,:) - plot_1_SEM(3,:));
%     x2 = [mvgc.parameters.freqs(1:mvgc.parameters.fres)', fliplr(mvgc.parameters.freqs(1:mvgc.parameters.fres)')];
%     inBetween = [curve1, fliplr(curve2)];
%     fill(x2, inBetween, 'k','FaceAlpha',0.6,'EdgeColor','none');
%     hold on
%     plot(mvgc.parameters.freqs(1:mvgc.parameters.fres),squeeze(plot_1_mean(3,:)),'linew', 2,'Color','k')
% 
%     curve1 = squeeze(plot_1_mean(4,:) + plot_1_SEM(4,:));
%     curve2 = squeeze(plot_1_mean(4,:) - plot_1_SEM(4,:));
%     x2 = [mvgc.parameters.freqs(1:mvgc.parameters.fres)', fliplr(mvgc.parameters.freqs(1:mvgc.parameters.fres)')];
%     inBetween = [curve1, fliplr(curve2)];
%     fill(x2, inBetween, 'k','FaceAlpha',0.2,'EdgeColor','none');
%     plot(mvgc.parameters.freqs(1:mvgc.parameters.fres),squeeze(plot_1_mean(4,:)),'linew', 2,'Color',[ 0 0 0 .5])
%     xline(6)
%     xlim([2 15])
%     ylim([0 .40])
%     xlabel('(Hz)','FontSize',9), ylabel('Granger Prediction','FontSize',11)
%     legend('','PL --> HPC','','HPC --> PL')
%     title('PL <--> HPC')
%     box off
% 
%     subplot(3,3,3)
%     curve1 = squeeze(plot_1_mean(5,:) + plot_1_SEM(5,:));
%     curve2 = squeeze(plot_1_mean(5,:) - plot_1_SEM(5,:));
%     x2 = [mvgc.parameters.freqs(1:mvgc.parameters.fres)', fliplr(mvgc.parameters.freqs(1:mvgc.parameters.fres)')];
%     inBetween = [curve1, fliplr(curve2)];
%     fill(x2, inBetween, 'k','FaceAlpha',0.6,'EdgeColor','none');
%     hold on
%     plot(mvgc.parameters.freqs(1:mvgc.parameters.fres),squeeze(plot_1_mean(5,:)),'linew', 2,'Color','k')
% 
%     curve1 = squeeze(plot_1_mean(6,:) + plot_1_SEM(6,:));
%     curve2 = squeeze(plot_1_mean(6,:) - plot_1_SEM(6,:));
%     x2 = [mvgc.parameters.freqs(1:mvgc.parameters.fres)', fliplr(mvgc.parameters.freqs(1:mvgc.parameters.fres)')];
%     inBetween = [curve1, fliplr(curve2)];
%     fill(x2, inBetween, 'k','FaceAlpha',0.2,'EdgeColor','none');
%     plot(mvgc.parameters.freqs(1:mvgc.parameters.fres),squeeze(plot_1_mean(6,:)),'linew', 2,'Color',[ 0 0 0 .5])
%     xline(6)
%     xlim([2 15])
%     ylim([0 .40])
%     xlabel('(Hz)','FontSize',9), ylabel('Granger Prediction','FontSize',11)
%     legend('','IL --> HPC','','HPC --> IL')
%     title('IL <--> HPC')
%     box off
% 
% 
% 
% 
%     % Freezing
%     subplot(3,3,4)
%     curve1 = squeeze(plot_2_mean(1,:) + plot_2_SEM(1,:));
%     curve2 = squeeze(plot_2_mean(1,:) - plot_2_SEM(1,:));
%     x2 = [mvgc.parameters.freqs(1:mvgc.parameters.fres)', fliplr(mvgc.parameters.freqs(1:mvgc.parameters.fres)')];
%     inBetween = [curve1, fliplr(curve2)];
%     fill(x2, inBetween, [.6, 0, 0],'FaceAlpha',0.2,'EdgeColor','none');
%     hold on
%     plot(mvgc.parameters.freqs(1:mvgc.parameters.fres),squeeze(plot_2_mean(1,:)),'linew', 2,'Color',[.6, 0, 0])
% 
%     curve1 = squeeze(plot_2_mean(2,:) + plot_2_SEM(2,:));
%     curve2 = squeeze(plot_2_mean(2,:) - plot_2_SEM(2,:));
%     x2 = [mvgc.parameters.freqs(1:mvgc.parameters.fres)', fliplr(mvgc.parameters.freqs(1:mvgc.parameters.fres)')];
%     inBetween = [curve1, fliplr(curve2)];
%     fill(x2, inBetween, [1, .7, 0],'FaceAlpha',0.2,'EdgeColor','none');
%     plot(mvgc.parameters.freqs(1:mvgc.parameters.fres),squeeze(plot_2_mean(2,:)),'linew', 2,'Color',[1, .7, 0])
%     xline(6)
%     xlim([2 15])
%     ylim([0 .40])
%     xlabel('(Hz)','FontSize',9), ylabel({'Freezing'; 'Granger Prediction'},'FontSize',11)
%     legend('','PL --> IL','','IL --> PL')
%     box off
% 
%     subplot(3,3,5)
%     curve1 = squeeze(plot_2_mean(3,:) + plot_2_SEM(3,:));
%     curve2 = squeeze(plot_2_mean(3,:) - plot_2_SEM(3,:));
%     x2 = [mvgc.parameters.freqs(1:mvgc.parameters.fres)', fliplr(mvgc.parameters.freqs(1:mvgc.parameters.fres)')];
%     inBetween = [curve1, fliplr(curve2)];
%     fill(x2, inBetween, [.6, 0, 0],'FaceAlpha',0.2,'EdgeColor','none');
%     hold on
%     plot(mvgc.parameters.freqs(1:mvgc.parameters.fres),squeeze(plot_2_mean(3,:)),'linew', 2,'Color',[.6, 0, 0])
% 
%     curve1 = squeeze(plot_2_mean(4,:) + plot_2_SEM(4,:));
%     curve2 = squeeze(plot_2_mean(4,:) - plot_2_SEM(4,:));
%     x2 = [mvgc.parameters.freqs(1:mvgc.parameters.fres)', fliplr(mvgc.parameters.freqs(1:mvgc.parameters.fres)')];
%     inBetween = [curve1, fliplr(curve2)];
%     fill(x2, inBetween, [1, .7, 0],'FaceAlpha',0.2,'EdgeColor','none');
%     plot(mvgc.parameters.freqs(1:mvgc.parameters.fres),squeeze(plot_2_mean(4,:)),'linew', 2,'Color',[1, .7, 0])
%     xline(6)
%     xlim([2 15])
%     ylim([0 .40])
%     xlabel('(Hz)','FontSize',9), ylabel('Granger Prediction','FontSize',11)
%     legend('','PL --> HPC','','HPC --> PL')
%     box off
% 
%     subplot(3,3,6)
%     curve1 = squeeze(plot_2_mean(5,:) + plot_2_SEM(5,:));
%     curve2 = squeeze(plot_2_mean(5,:) - plot_2_SEM(5,:));
%     x2 = [mvgc.parameters.freqs(1:mvgc.parameters.fres)', fliplr(mvgc.parameters.freqs(1:mvgc.parameters.fres)')];
%     inBetween = [curve1, fliplr(curve2)];
%     fill(x2, inBetween, [.6, 0, 0],'FaceAlpha',0.2,'EdgeColor','none');
%     hold on
%     plot(mvgc.parameters.freqs(1:mvgc.parameters.fres),squeeze(plot_2_mean(5,:)),'linew', 2,'Color',[.6, 0, 0])
% 
%     curve1 = squeeze(plot_2_mean(6,:) + plot_2_SEM(6,:));
%     curve2 = squeeze(plot_2_mean(6,:) - plot_2_SEM(6,:));
%     x2 = [mvgc.parameters.freqs(1:mvgc.parameters.fres)', fliplr(mvgc.parameters.freqs(1:mvgc.parameters.fres)')];
%     inBetween = [curve1, fliplr(curve2)];
%     fill(x2, inBetween, [1, .7, 0],'FaceAlpha',0.2,'EdgeColor','none');
%     plot(mvgc.parameters.freqs(1:mvgc.parameters.fres),squeeze(plot_2_mean(6,:)),'linew', 2,'Color',[1, .7, 0])
%     xline(6)
%     xlim([2 15])
%     ylim([0 .40])
%     xlabel('(Hz)','FontSize',9), ylabel('Granger Prediction','FontSize',11)
%     legend('','IL --> HPC','','HPC --> IL')
%     box off
% 
% 
% 
%     
%     % Non-Freezing
%     subplot(3,3,7)
%     curve1 = squeeze(plot_3_mean(1,:) + plot_3_SEM(1,:));
%     curve2 = squeeze(plot_3_mean(1,:) - plot_3_SEM(1,:));
%     x2 = [mvgc.parameters.freqs(1:mvgc.parameters.fres)', fliplr(mvgc.parameters.freqs(1:mvgc.parameters.fres)')];
%     inBetween = [curve1, fliplr(curve2)];
%     fill(x2, inBetween, [.2, .2, 1],'FaceAlpha',0.2,'EdgeColor','none');
%     hold on
%     plot(mvgc.parameters.freqs(1:mvgc.parameters.fres),squeeze(plot_3_mean(1,:)),'linew', 2,'Color',[.2, .2, 1])
% 
%     curve1 = squeeze(plot_3_mean(2,:) + plot_3_SEM(2,:));
%     curve2 = squeeze(plot_3_mean(2,:) - plot_3_SEM(2,:));
%     x2 = [mvgc.parameters.freqs(1:mvgc.parameters.fres)', fliplr(mvgc.parameters.freqs(1:mvgc.parameters.fres)')];
%     inBetween = [curve1, fliplr(curve2)];
%     fill(x2, inBetween, [.2, .8, 1],'FaceAlpha',0.2,'EdgeColor','none');
%     plot(mvgc.parameters.freqs(1:mvgc.parameters.fres),squeeze(plot_3_mean(2,:)),'linew', 2,'Color',[.2, .8, 1])
%     xline(6)
%     xlim([2 15])
%     ylim([0 .40])
%     xlabel('(Hz)','FontSize',9), ylabel({'Movement activity'; 'Granger Prediction'},'FontSize',11)
%     legend('','PL --> IL','','IL --> PL')
%     box off
% 
%     subplot(3,3,8)
%     curve1 = squeeze(plot_3_mean(3,:) + plot_3_SEM(3,:));
%     curve2 = squeeze(plot_3_mean(3,:) - plot_3_SEM(3,:));
%     x2 = [mvgc.parameters.freqs(1:mvgc.parameters.fres)', fliplr(mvgc.parameters.freqs(1:mvgc.parameters.fres)')];
%     inBetween = [curve1, fliplr(curve2)];
%     fill(x2, inBetween, [.2, .2, 1],'FaceAlpha',0.2,'EdgeColor','none');
%     hold on
%     plot(mvgc.parameters.freqs(1:mvgc.parameters.fres),squeeze(plot_3_mean(3,:)),'linew', 2,'Color',[.2, .2, 1])
% 
%     curve1 = squeeze(plot_3_mean(4,:) + plot_3_SEM(4,:));
%     curve2 = squeeze(plot_3_mean(4,:) - plot_3_SEM(4,:));
%     x2 = [mvgc.parameters.freqs(1:mvgc.parameters.fres)', fliplr(mvgc.parameters.freqs(1:mvgc.parameters.fres)')];
%     inBetween = [curve1, fliplr(curve2)];
%     fill(x2, inBetween, [.2, .8, 1],'FaceAlpha',0.2,'EdgeColor','none');
%     plot(mvgc.parameters.freqs(1:mvgc.parameters.fres),squeeze(plot_3_mean(4,:)),'linew', 2,'Color',[.2, .8, 1])
%     xline(6)
%     xlim([2 15])
%     ylim([0 .40])
%     xlabel('(Hz)','FontSize',9), ylabel('Granger Prediction','FontSize',11)
%     legend('','PL --> HPC','','HPC --> PL')
%     box off
% 
%     subplot(3,3,9)
%     curve1 = squeeze(plot_3_mean(5,:) + plot_3_SEM(5,:));
%     curve2 = squeeze(plot_3_mean(5,:) - plot_3_SEM(5,:));
%     x2 = [mvgc.parameters.freqs(1:mvgc.parameters.fres)', fliplr(mvgc.parameters.freqs(1:mvgc.parameters.fres)')];
%     inBetween = [curve1, fliplr(curve2)];
%     fill(x2, inBetween, [.2, .2, 1],'FaceAlpha',0.2,'EdgeColor','none');
%     hold on
%     plot(mvgc.parameters.freqs(1:mvgc.parameters.fres),squeeze(plot_3_mean(5,:)),'linew', 2,'Color',[.2, .2, 1])
% 
%     curve1 = squeeze(plot_3_mean(6,:) + plot_3_SEM(6,:));
%     curve2 = squeeze(plot_3_mean(6,:) - plot_3_SEM(6,:));
%     x2 = [mvgc.parameters.freqs(1:mvgc.parameters.fres)', fliplr(mvgc.parameters.freqs(1:mvgc.parameters.fres)')];
%     inBetween = [curve1, fliplr(curve2)];
%     fill(x2, inBetween, [.2, .8, 1],'FaceAlpha',0.2,'EdgeColor','none');
%     plot(mvgc.parameters.freqs(1:mvgc.parameters.fres),squeeze(plot_3_mean(6,:)),'linew', 2,'Color',[.2, .8, 1])
%     xline(6)
%     xlim([2 15])
%     ylim([0 .40])
%     xlabel('(Hz)','FontSize',9), ylabel('Granger Prediction','FontSize',11)
%     legend('','IL --> HPC','','HPC --> IL')
%     box off
% 
% end
% 
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
