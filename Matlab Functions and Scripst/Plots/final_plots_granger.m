% Final Plots Granger

% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  01/2024
% Last update: 02/2024

% Extinction - CS-Trials without noise
% CSIT_extinction{1} = 1:5;      % MT1
% CSIT_extinction{2} = 1:4;      % MT3
% CSIT_extinction{3} = [1 2 4];  % MT4
% CSIT_extinction{4} = 1:5;      % MT5
% CSIT_extinction{5} = 1:5;      % MT6
% CSIT_extinction{6} = 1:5;      % MT7

% % Retrieval - CS-Trials without noise
% CSIT_retrieval{1} = 1:5;         % MT1
% CSIT_retrieval{2} = [1 3 4 5];   % MT3
% CSIT_retrieval{3} = [1];         % MT4
% CSIT_retrieval{4} = 1:5;         % MT5
% CSIT_retrieval{5} = 1:5;         % MT6
% CSIT_retrieval{6} = 1:5;         % MT7



%% Prepara data

%% Plot specs version 3
%Lazy code... I'm dead...


%% Loading data - CS
%data_2_plot_CS_extinction(:,:,:,ms) = mvgc.stats.Spect_2_12Hz_2_plot_mean;
%data_2_plot_CS_retrieval(:,:,:,ms) = mvgc.stats.Spect_2_12Hz_2_plot_mean;


%% Loading data - Baseline and Post-Tone
% data_2_plot_extinction_baseline(:,:,ms) = mvgc.stats.Spect_2_12Hz_baseline_2_plot;
% data_2_plot_extinction_post_tone(:,:,ms) = mvgc.stats.Spect_2_12Hz_post_tone_2_plot;

data_2_plot_retrieval_baseline(:,:,ms) = mvgc.stats.Spect_2_12Hz_baseline_2_plot;
data_2_plot_retrieval_post_tone(:,:,ms) = mvgc.stats.Spect_2_12Hz_post_tone_2_plot;

%% CS - Mean and SEM

% Extinction - Baseline - CS
% data_2_plot_CS_extinction_mean = mean(data_2_plot_CS_extinction,4);
% data_2_plot_CS_extinction_SEM  = std(data_2_plot_CS_extinction,[],4)./size(data_2_plot_CS_extinction,4);

%retrieval - Baseline - CS
% data_2_plot_CS_retrieval_mean  = mean(data_2_plot_CS_retrieval,4);
% data_2_plot_CS_retrieval_SEM   = std(data_2_plot_CS_retrieval,[],4)./size(data_2_plot_CS_retrieval,4);

% freq_v = mvgc.parameters.freqs(mvgc.parameters.frex_idx_2_12Hz);

%% Baseline pos-tone - - Mean and SEM

% Extinction Baseline
data_2_plot_extinction_baseline_mean = mean(data_2_plot_extinction_baseline(:,:,[1 2 3 4 5 6]),3);
data_2_plot_extinction_baseline_SEM  = std(data_2_plot_extinction_baseline(:,:,[1 2 3 4 5 6]),[],3)./size(data_2_plot_extinction_baseline(:,:,[1 3 4 5 6]),3);

% Extinction Post-Tone
data_2_plot_extinction_post_tone_mean = mean(data_2_plot_extinction_post_tone(:,:,[1 2 3 4 5 6]),3);
data_2_plot_extinction_post_tone_SEM  = std(data_2_plot_extinction_post_tone(:,:,[1 2 3 4 5 6]),[],3)./size(data_2_plot_extinction_post_tone(:,:,[1 3 4 5 6]),3);

% % Retrieval
data_2_plot_retrieval_baseline_mean = mean(data_2_plot_retrieval_baseline(:,:,[1 2 3 4 5 6]),3);
data_2_plot_retrieval_baseline_SEM  = std(data_2_plot_retrieval_baseline(:,:,[1 2 3 4 5 6]),[],3)./size(data_2_plot_retrieval_baseline(:,:,[1 2 3 4 5 6]),3);

% % Extinction Post-Tone
data_2_plot_retrieval_post_tone_mean = mean(data_2_plot_retrieval_post_tone(:,:,[1 2 3 4 5 6]),3);
data_2_plot_retrieval_post_tone_SEM  = std(data_2_plot_retrieval_post_tone(:,:,[1 2 3 4 5 6]),[],3)./size(data_2_plot_retrieval_post_tone(:,:,[1 2 3 4 5 6]),3);

% 
 freq_v = mvgc.parameters.freqs(mvgc.parameters.frex_idx_2_12Hz);
%% Choose data toplot 

% CS
% plot_1_mean = squeeze(data_2_plot_extinction_mean(:,2,:));
% plot_1_SEM  = squeeze(data_2_plot_extinction_SEM(:,2,:));
% 
% plot_2_mean = squeeze(data_2_plot_retrieval_mean(:,2,:));
% plot_2_SEM  = squeeze(data_2_plot_retrieval_SEM(:,2,:));

%%
% Baseline and Post-Tone
plot_1_mean = data_2_plot_extinction_baseline_mean;
plot_1_SEM  = data_2_plot_extinction_baseline_SEM;

plot_2_mean = data_2_plot_extinction_post_tone_mean;
plot_2_SEM  = data_2_plot_extinction_post_tone_SEM;

plot_1_mean = data_2_plot_retrieval_baseline_mean;
plot_1_SEM  = data_2_plot_retrieval_baseline_SEM;

plot_2_mean = data_2_plot_retrieval_post_tone_mean;
plot_2_SEM  = data_2_plot_retrieval_post_tone_SEM;

%%
figure
set(gcf,'color','w');
set(gcf, 'Position', get(0, 'Screensize'));

sgtitle({['Multivariate Granger Causality - VAR model estimation regression mode (LWR). Order ~= 35'];[]})

% CS-Trials

subplot(2,3,1)

% SEM shades
curve1 = squeeze(plot_1_mean(1,:) + plot_1_SEM(1,:));
curve2 = squeeze(plot_1_mean(1,:) - plot_1_SEM(1,:));
x2 = [mvgc.parameters.freqs', fliplr(mvgc.parameters.freqs')];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [.6, 0, 0],'FaceAlpha',0.2,'EdgeColor','none');
hold on
plot(mvgc.parameters.freqs,squeeze(plot_1_mean(1,:)),'linew', 2,'Color',[.6, 0, 0])

curve1 = squeeze(plot_1_mean(2,:) + plot_1_SEM(2,:));
curve2 = squeeze(plot_1_mean(2,:) - plot_1_SEM(2,:));
x2 = [mvgc.parameters.freqs', fliplr(mvgc.parameters.freqs')];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [1, .7, 0],'FaceAlpha',0.2,'EdgeColor','none');
plot(mvgc.parameters.freqs,squeeze(plot_1_mean(2,:)),'linew', 2,'Color',[1, .7, 0])

xline(6)
xlim([2 15])
ylim([0 .3])
xlabel('(Hz)','FontSize',9), ylabel('Granger Prediction','FontSize',11)
legend('','PL --> IL','','IL --> PL')
title('PL <--> IL')
box off

subplot(2,3,2)
curve1 = squeeze(plot_1_mean(3,:) + plot_1_SEM(3,:));
curve2 = squeeze(plot_1_mean(3,:) - plot_1_SEM(3,:));
x2 = [mvgc.parameters.freqs', fliplr(mvgc.parameters.freqs')];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [.6, 0, 0],'FaceAlpha',0.2,'EdgeColor','none');
hold on
plot(mvgc.parameters.freqs,squeeze(plot_1_mean(3,:)),'linew', 2,'Color',[.6, 0, 0])

curve1 = squeeze(plot_1_mean(4,:) + plot_1_SEM(4,:));
curve2 = squeeze(plot_1_mean(4,:) - plot_1_SEM(4,:));
x2 = [mvgc.parameters.freqs', fliplr(mvgc.parameters.freqs')];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [1, .7, 0],'FaceAlpha',0.2,'EdgeColor','none');
plot(mvgc.parameters.freqs,squeeze(plot_1_mean(4,:)),'linew', 2,'Color',[1, .7, 0])
xline(6)
xlim([2 15])
ylim([0 .3])
xlabel('(Hz)','FontSize',9), ylabel('Granger Prediction','FontSize',11)
legend('','PL --> HPC','','HPC --> PL')
title('PL <--> HPC')
box off

subplot(2,3,3)
curve1 = squeeze(plot_1_mean(5,:) + plot_1_SEM(5,:));
curve2 = squeeze(plot_1_mean(5,:) - plot_1_SEM(5,:));
x2 = [mvgc.parameters.freqs', fliplr(mvgc.parameters.freqs')];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [.6, 0, 0],'FaceAlpha',0.2,'EdgeColor','none');
hold on
plot(mvgc.parameters.freqs,squeeze(plot_1_mean(5,:)),'linew', 2,'Color',[.6, 0, 0])

curve1 = squeeze(plot_1_mean(6,:) + plot_1_SEM(6,:));
curve2 = squeeze(plot_1_mean(6,:) - plot_1_SEM(6,:));
x2 = [mvgc.parameters.freqs', fliplr(mvgc.parameters.freqs')];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [1, .7, 0],'FaceAlpha',0.2,'EdgeColor','none');
plot(mvgc.parameters.freqs,squeeze(plot_1_mean(6,:)),'linew', 2,'Color',[1, .7, 0])
xline(6)
xlim([2 15])
ylim([0 .3])
xlabel('(Hz)','FontSize',9), ylabel('Granger Prediction','FontSize',11)
legend('','IL --> HPC','','HPC --> IL')
title('IL <--> HPC')
box off

subplot(2,3,4)
curve1 = squeeze(plot_2_mean(1,:) + plot_2_SEM(1,:));
curve2 = squeeze(plot_2_mean(1,:) - plot_2_SEM(1,:));
x2 = [mvgc.parameters.freqs', fliplr(mvgc.parameters.freqs')];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [.2, .2, 1],'FaceAlpha',0.2,'EdgeColor','none');
hold on
plot(mvgc.parameters.freqs,squeeze(plot_2_mean(1,:)),'linew', 2,'Color',[.2, .2, 1])

curve1 = squeeze(plot_2_mean(2,:) + plot_2_SEM(2,:));
curve2 = squeeze(plot_2_mean(2,:) - plot_2_SEM(2,:));
x2 = [mvgc.parameters.freqs', fliplr(mvgc.parameters.freqs')];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [.2, .8, 1],'FaceAlpha',0.2,'EdgeColor','none');
plot(mvgc.parameters.freqs,squeeze(plot_2_mean(2,:)),'linew', 2,'Color',[.2, .8, 1])
xline(6)
xlim([2 15])
ylim([0 .3])
xlabel('(Hz)','FontSize',9), ylabel('Granger Prediction','FontSize',11)
legend('','PL --> IL','','IL --> PL')
box off

subplot(2,3,5)
curve1 = squeeze(plot_2_mean(3,:) + plot_2_SEM(3,:));
curve2 = squeeze(plot_2_mean(3,:) - plot_2_SEM(3,:));
x2 = [mvgc.parameters.freqs', fliplr(mvgc.parameters.freqs')];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [.2, .2, 1],'FaceAlpha',0.2,'EdgeColor','none');
hold on
plot(mvgc.parameters.freqs,squeeze(plot_2_mean(3,:)),'linew', 2,'Color',[.2, .2, 1])

curve1 = squeeze(plot_2_mean(4,:) + plot_2_SEM(4,:));
curve2 = squeeze(plot_2_mean(4,:) - plot_2_SEM(4,:));
x2 = [mvgc.parameters.freqs', fliplr(mvgc.parameters.freqs')];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [.2, .8, 1],'FaceAlpha',0.2,'EdgeColor','none');
plot(mvgc.parameters.freqs,squeeze(plot_2_mean(4,:)),'linew', 2,'Color',[.2, .8, 1])
xline(6)
xlim([2 15])
ylim([0 .3])
xlabel('(Hz)','FontSize',9), ylabel('Granger Prediction','FontSize',11)
legend('','PL --> HPC','','HPC --> PL')
box off

subplot(2,3,6)
curve1 = squeeze(plot_2_mean(5,:) + plot_2_SEM(5,:));
curve2 = squeeze(plot_2_mean(5,:) - plot_2_SEM(5,:));
x2 = [mvgc.parameters.freqs', fliplr(mvgc.parameters.freqs')];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [.2, .2, 1],'FaceAlpha',0.2,'EdgeColor','none');
hold on
plot(mvgc.parameters.freqs,squeeze(plot_2_mean(5,:)),'linew', 2,'Color',[.2, .2, 1])

curve1 = squeeze(plot_2_mean(6,:) + plot_2_SEM(6,:));
curve2 = squeeze(plot_2_mean(6,:) - plot_2_SEM(6,:));
x2 = [mvgc.parameters.freqs', fliplr(mvgc.parameters.freqs')];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [.2, .8, 1],'FaceAlpha',0.2,'EdgeColor','none');
plot(mvgc.parameters.freqs,squeeze(plot_2_mean(6,:)),'linew', 2,'Color',[.2, .8, 1])
xline(6)
xlim([2 15])
ylim([0 .3])
xlabel('(Hz)','FontSize',9), ylabel('Granger Prediction','FontSize',11)
legend('','IL --> HPC','','HPC --> IL')
box off

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
% 
% 
% %% plot Matrix - 3 - 6 Hz
% % Plot Retrival session at the second line
% 
% data_2_plot_matrix{1,1} = mean(cell2mat(mvgc.Fint_3_6Hz(1,1)),3); % baseline
% data_2_plot_matrix{1,2} = mean(reshape(cell2mat(mvgc.Fint_3_6Hz(2,CSIT{1})),size(mvgc.Fint_3_6Hz,1),[],length(CSIT{1})),3); % CS-Trial
% data_2_plot_matrix{1,3} = mean(reshape(cell2mat(mvgc.Fint_3_6Hz(3,CSIT{1})),size(mvgc.Fint_3_6Hz,1),[],length(CSIT{1})),3); % ITI
% 
% data_2_plot_matrix{2,1} = mean(cell2mat(mvgc.Fint_6_9Hz(1,1)),3); % baseline
% data_2_plot_matrix{2,2} = mean(reshape(cell2mat(mvgc.Fint_6_9Hz(2,CSIT{1})),size(mvgc.Fint_6_9Hz,1),[],length(CSIT{1})),3); % CS-Trial
% data_2_plot_matrix{2,3} = mean(reshape(cell2mat(mvgc.Fint_6_9Hz(3,CSIT{1})),size(mvgc.Fint_6_9Hz,1),[],length(CSIT{1})),3); % ITI
% 
% %%
% figure
% set(gcf,'color','w');
% set(gcf, 'Position', get(0, 'Screensize'));
% 
% myLabel = {'mPFC PL','mPFC IL','HPC'};
% 
% sgtitle( '3-6 Hz')
% 
% subplot(231)
% imagesc(data_2_plot_matrix{1,1},'AlphaData',~isnan(data_2_plot_matrix{1,1}))
% colorbar
% colormap jet
% xticklabels(myLabel)
% xticks(1:length(myLabel))
% %xtickangle(90)
% yticks(1:length(myLabel))
% yticklabels(myLabel)
% xlabel('- From - ','FontSize',14), ylabel('- To - ','FontSize',14)
% title('Baseline','FontSize',12)
% caxis([0 .15])
% 
% subplot(232)
% imagesc(data_2_plot_matrix{1,2},'AlphaData',~isnan(data_2_plot_matrix{1,2}))
% colorbar
% colormap jet
% xticklabels(myLabel)
% xticks(1:length(myLabel))
% %xtickangle(90)
% yticks(1:length(myLabel))
% yticklabels(myLabel)
% xlabel('- From - ','FontSize',14), ylabel('- To - ','FontSize',14)
% title('CS-Trials','FontSize',12)
% caxis([0 .15])
% 
% subplot(233)
% imagesc(data_2_plot_matrix{1,3},'AlphaData',~isnan(data_2_plot_matrix{1,3}))
% colorbar
% colormap jet
% xticklabels(myLabel)
% xticks(1:length(myLabel))
% %xtickangle(90)
% yticks(1:length(myLabel))
% yticklabels(myLabel)
% xlabel('- From - ','FontSize',14), ylabel('- To - ','FontSize',14)
% title('ITI','FontSize',12)
% caxis([0 .15])
% 
% %% Save
% newStr1 = id(1:end-8);
% name_1  = strcat(Path,'/',newStr1,'_Granger_Spec_Mean_All_3');
% Path    = files.FilesLoaded{1,1}(ms).folder;
% saveas(gcf,name_1,'png') % save figure
% 
% close all
% 
% clear('ii','jj')
% clear('Combinations_2','Combinations_1','sub_idx','steps','data2plot_mean','cc','session','name_1','newStr1','path','data2plot') 
% 
% 
% %% plot Matrix - 6 - 9 Hz
% % Plot Retrival session at the second line
% 
% figure
% set(gcf,'color','w');
% set(gcf, 'Position', get(0, 'Screensize'));
% 
% sgtitle( '6-9 Hz')
% subplot(231)
% imagesc(data_2_plot_matrix{2,1},'AlphaData',~isnan(data_2_plot_matrix{1,1}))
% colorbar
% colormap jet
% xticklabels(myLabel)
% xticks(1:length(myLabel))
% %xtickangle(90)
% yticks(1:length(myLabel))
% yticklabels(myLabel)
% xlabel('- From - ','FontSize',14), ylabel('- To - ','FontSize',14)
% title('Baseline','FontSize',12)
% caxis([0 .15])
% 
% subplot(232)
% imagesc(data_2_plot_matrix{2,2},'AlphaData',~isnan(data_2_plot_matrix{1,2}))
% colorbar
% colormap jet
% xticklabels(myLabel)
% xticks(1:length(myLabel))
% %xtickangle(90)
% yticks(1:length(myLabel))
% yticklabels(myLabel)
% xlabel('- From - ','FontSize',14), ylabel('- To - ','FontSize',14)
% title('CS-Trials','FontSize',12)
% caxis([0 .15])
% 
% subplot(233)
% imagesc(data_2_plot_matrix{2,3},'AlphaData',~isnan(data_2_plot_matrix{1,3}))
% colorbar
% colormap jet
% xticklabels(myLabel)
% xticks(1:length(myLabel))
% %xtickangle(90)
% yticks(1:length(myLabel))
% yticklabels(myLabel)
% xlabel('- From - ','FontSize',14), ylabel('- To - ','FontSize',14)
% title('ITI','FontSize',12)
% caxis([0 .15])
% 
% %% Save
% newStr1 = id(1:end-8);
% name_1  = strcat(Path,'/',newStr1,'_Granger_Spec_Mean_All_3');
% Path    = files.FilesLoaded{1,1}(ms).folder;
% saveas(gcf,name_1,'png') % save figure
% 
% close all
% 
% clear('ii','jj')
% clear('Combinations_2','Combinations_1','sub_idx','steps','data2plot_mean','cc','session','name_1','newStr1','path','data2plot') 
