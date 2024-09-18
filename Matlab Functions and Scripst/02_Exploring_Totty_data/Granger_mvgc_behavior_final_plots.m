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
% 
for ii = 1:size(mvgc.stats.spectrum_2_12Hz_mean,1)
    mvgc.stats.spectrum_2_12Hz_mean{ii,2}(isempty(mvgc.stats.spectrum_2_12Hz_mean{ii,2})) = NaN;
    data_2_plot_extinction{ii,1}(:,:,:,ms) = mvgc.stats.spectrum_2_12Hz_mean{ii,2};
end

freq_v = mvgc.parameters.freqs(mvgc.parameters.frex_idx_2_12Hz);


%% Load data Retrieval

% data organization from mvgc.F_spect --> see data_preprocessing.m
% - Row 1: Baseline freezing
% - Row 2: Baseline not freezing
% - Row 3: CS-TONE freezing
% - Row 4: CS-TONE not freezing
% - Row 5: ITI freezing
% - Row 6: ITI not freezing

% for ii = 1:size(mvgc.stats.spectrum_2_12Hz_mean,1)
%     mvgc.stats.spectrum_2_12Hz_mean{ii,2}(isempty(mvgc.stats.spectrum_2_12Hz_mean{ii,2})) = NaN;
%     data_2_plot_retrieval{ii,1}(:,:,:,ms) = mvgc.stats.spectrum_2_12Hz_mean{ii,2};
% end
% 
% % Frequency vector
% 
% freq_v = mvgc.parameters.freqs(mvgc.parameters.frex_idx_2_12Hz);

%%

%% Mean and SEM

% 
% data_2_plot_extinction_mean_SEM = [];
% 
% for ii = 1:size(data_2_plot_extinction,1)
%     data_2_plot_extinction_mean_SEM{ii,1}  = 10.*mean(data_2_plot_extinction{ii,1},4,'omitnan');
%     data_2_plot_extinction_mean_SEM{ii,2}  = 10.*std(data_2_plot_extinction{ii,1},[],4,'omitnan')./size(data_2_plot_extinction{ii,1},4);
%     
%     % SEM shades
%     curve1_ext{ii,1} = data_2_plot_extinction_mean_SEM{ii, 1} + data_2_plot_extinction_mean_SEM{ii, 2};
%     curve2_ext{ii,1} = data_2_plot_extinction_mean_SEM{ii, 1} - data_2_plot_extinction_mean_SEM{ii, 2};
%     data_2_plot_extinction_mean_SEM{ii,3} = cat(3,curve1_ext{ii,1},flip(curve2_ext{ii,1},3));
%     data_2_plot_extinction_mean_SEM{ii,4} = [freq_v', fliplr(freq_v')];
% 
% end
% 
% % 
% % data_2_plot_retrieval_mean_SEM = [];
% % 
% % for ii = 1:size(data_2_plot_retrieval,1)
% %     data_2_plot_retrieval_mean_SEM{ii,1}  = 10.*mean(data_2_plot_retrieval{ii,1},4,'omitnan');
% %     data_2_plot_retrieval_mean_SEM{ii,2}  = 10.*std(data_2_plot_retrieval{ii,1},[],4,'omitnan')./size(data_2_plot_retrieval{ii,1},4);
%     
%     % SEM shades
%     curve1_ret{ii,1} = data_2_plot_retrieval_mean_SEM{ii, 1} + data_2_plot_retrieval_mean_SEM{ii, 2};
%     curve2_ret{ii,1} = data_2_plot_retrieval_mean_SEM{ii, 1} - data_2_plot_retrieval_mean_SEM{ii, 2};
%     data_2_plot_retrieval_mean_SEM{ii,3} = cat(3,curve1_ret{ii,1},flip(curve2_ret{ii,1},3));
%     data_2_plot_retrieval_mean_SEM{ii,4} = [freq_v', fliplr(freq_v')];
% 
% end


%% Plot
% Plot spectral causal graph. Version 2.
% Plot pairwise spectral quantities in |P|, a 3-dim numerical matrix with
% first index representing target ("to"), second index source ("from")
% and third index frequency range - typically spectral causalities

% all possible combinations
Combinations_2 = nchoosek(1:3,2);
Combinations_1 = flip(Combinations_2,2); % all possible combinations

% data to Plot
data_2_plot = data_2_plot_extinction_mean_SEM;
%data_2_plot = data_2_plot_retrieval_mean_SEM;

figure
set(gcf,'color','w');
sc = [1,1,960,1200];
set(gcf, 'Position', sc);
sgtitle({'Pairwise-conditional Granger causality - frequency domain.';''});


% Baseline
for cc = 1:size(Combinations_1,1)

    subplot(6,3,cc)
    fill(data_2_plot{1,4},squeeze(data_2_plot{1,3}(Combinations_1(cc,1),Combinations_1(cc,2),:)), 'k','FaceAlpha',0.5,'EdgeColor','none')
    hold on
    plot(freq_v,squeeze(data_2_plot{1,1}(Combinations_1(cc,1),Combinations_1(cc,2),:)),'k','linew',2)
    fill(data_2_plot{1,4},squeeze(data_2_plot{1,3}(Combinations_2(cc,1),Combinations_2(cc,2),:)), 'k','FaceAlpha',0.2,'EdgeColor','none')
    plot(freq_v,squeeze(data_2_plot{1,1}(Combinations_2(cc,1),Combinations_2(cc,2),:)),'k','linew',2)
   

    if          cc == 1
        legend('PL --> IL','','IL --> PL','')
        legend('boxoff')
        legend('FontSize',8)
        ylabel({'Baseline - Freezing';'Granger causality'},'FontWeight','bold','FontSize',12)

    elseif      cc == 2
        legend('PL --> HPC','','HPC --> PL','')
        legend('boxoff')
        legend('FontSize',8)

    else
        legend('IL --> HPC','','HPC --> IL','')
        legend('boxoff')
        legend('FontSize',8)

    end

    xlabel('(Hz)','FontSize',12)
    xlim([2 11])
    ylim([0.005 6*10^-2])

end



for cc = 1:size(Combinations_1,1)

    subplot(6,3,cc+3)
    fill(data_2_plot{2,4},squeeze(data_2_plot{2,3}(Combinations_1(cc,1),Combinations_1(cc,2),:)),'k','FaceAlpha',0.3,'EdgeColor','none')
    hold on
    plot(freq_v,squeeze(data_2_plot{2,1}(Combinations_1(cc,1),Combinations_1(cc,2),:)),'k','linew',2)
    fill(data_2_plot{2,4},squeeze(data_2_plot{2,3}(Combinations_2(cc,1),Combinations_2(cc,2),:)),[.5 .5 .5],'FaceAlpha',0.2,'EdgeColor','none')
    plot(freq_v,squeeze(data_2_plot{2,1}(Combinations_2(cc,1),Combinations_2(cc,2),:)),'color',[.5 .5 .5],'linew',2)

    if          cc == 1
        legend('PL --> IL','', 'IL --> PL','')
        legend('boxoff')
        legend('FontSize',8)       
        ylabel({'Baseline - NonFreezing';'Granger causality'},'FontWeight','bold','FontSize',12)
    elseif      cc == 2
        legend('PL --> HPC','','HPC --> PL','')
        legend('boxoff')
        legend('FontSize',8)        
    else
        legend('IL --> HPC','', 'HPC --> IL','')
        legend('boxoff')
        legend('FontSize',8)        
    end

    xlabel('(Hz)','FontSize',12)
    xlim([2 11])
    ylim([0.005 6*10^-2])

end



%CS-Tone
for cc = 1:size(Combinations_1,1)

    subplot(6,3,cc+6)
    fill(data_2_plot{3,4},squeeze(data_2_plot{3,3}(Combinations_1(cc,1),Combinations_1(cc,2),:)),[.6, 0, 0],'FaceAlpha',0.5,'EdgeColor','none')
    hold on
    plot(freq_v,squeeze(data_2_plot{3,1}(Combinations_1(cc,1),Combinations_1(cc,2),:)),'color',[.6, 0, 0],'linew',2)
    fill(data_2_plot{3,4},squeeze(data_2_plot{3,3}(Combinations_2(cc,1),Combinations_2(cc,2),:)),[1, 0, .4],'FaceAlpha',0.2,'EdgeColor','none')
    plot(freq_v,squeeze(data_2_plot{3,1}(Combinations_2(cc,1),Combinations_2(cc,2),:)),'color',[1, 0, .4],'linew',2)

    if          cc == 1
        legend('PL --> IL','', 'IL --> PL','')
        legend('boxoff')
        legend('FontSize',8)        
        ylabel({'CS-Tone Freezing';'Granger causality'},'FontWeight','bold','FontSize',12)
    elseif      cc == 2
        legend('PL --> HPC','','HPC --> PL','')
        legend('boxoff')
        legend('FontSize',8)        
    else
        legend('IL --> HPC','', 'HPC --> IL','')
        legend('boxoff')
        legend('FontSize',8)        
    end

    xlabel('(Hz)','FontSize',12)
    xlim([2 11])
    ylim([0.005 6*10^-2])

end

for cc = 1:size(Combinations_1,1)

    subplot(6,3,cc+9)
    fill(data_2_plot{4,4},squeeze(data_2_plot{4,3}(Combinations_1(cc,1),Combinations_1(cc,2),:)),[1, .8, 0],'FaceAlpha',0.5,'EdgeColor','none')
    hold on
    plot(freq_v,squeeze(data_2_plot{4,1}(Combinations_1(cc,1),Combinations_1(cc,2),:)),'Color',[1, .4, 0],'linew',2)
    fill(data_2_plot{4,4},squeeze(data_2_plot{4,3}(Combinations_2(cc,1),Combinations_2(cc,2),:)),[1, .8, .4],'FaceAlpha',0.2,'EdgeColor','none')
    plot(freq_v,squeeze(data_2_plot{4,1}(Combinations_2(cc,1),Combinations_2(cc,2),:)),'Color',[1, .7, .0],'linew',2)

    if          cc == 1
        legend('PL --> IL','', 'IL --> PL','')
        legend('boxoff')
        legend('FontSize',8)        
        ylabel({'CS-Tone NonFreezing';'Granger causality'},'FontWeight','bold','FontSize',12)
    elseif      cc == 2
        legend('PL --> HPC','','HPC --> PL','')
        legend('boxoff')
        legend('FontSize',8)       
    else
        legend('IL --> HPC','', 'HPC --> IL','')
        legend('boxoff')
        legend('FontSize',8)

    end

    xlabel('(Hz)','FontSize',12)
    xlim([2 11])
    ylim([0.005 6*10^-2])

end


%ITI

for cc = 1:size(Combinations_1,1)

    subplot(6,3,cc+12)
    fill(data_2_plot{5,4},squeeze(data_2_plot{5,3}(Combinations_1(cc,1),Combinations_1(cc,2),:)),[.2, .2, 1],'FaceAlpha',0.5,'EdgeColor','none')
    hold on
    plot(freq_v,squeeze(data_2_plot{5,1}(Combinations_1(cc,1),Combinations_1(cc,2),:)),'Color',[.2, .2, 1],'linew',2)
    fill(data_2_plot{5,4},squeeze(data_2_plot{5,3}(Combinations_2(cc,1),Combinations_2(cc,2),:)),[.2, .2, 1],'FaceAlpha',0.2,'EdgeColor','none')
    plot(freq_v,squeeze(data_2_plot{5,1}(Combinations_2(cc,1),Combinations_2(cc,2),:)),'Color',[.2, .2, 1],'linew',2)

    if          cc == 1
        legend('PL --> IL','', 'IL --> PL','')
        legend('boxoff')
        legend('FontSize',8)       
        ylabel({'ITI Freezing';'Granger causality'},'FontWeight','bold','FontSize',12)
    elseif      cc == 2
        legend('PL --> HPC','','HPC --> PL','')
        legend('boxoff')
        legend('FontSize',8)   

    else
        legend('IL --> HPC','', 'HPC --> IL','')
        legend('boxoff')
        legend('FontSize',8)       
    end

    xlabel('(Hz)','FontSize',12)
    xlim([2 11])
    ylim([0.005 6*10^-2])

end

for cc = 1:size(Combinations_1,1)

    subplot(6,3,cc+15)
    fill(data_2_plot{6,4},squeeze(data_2_plot{6,3}(Combinations_1(cc,1),Combinations_1(cc,2),:)),[.2, .8, 1],'FaceAlpha',0.5,'EdgeColor','none')
    hold on
    plot(freq_v,squeeze(data_2_plot{6,1}(Combinations_1(cc,1),Combinations_1(cc,2),:)),'Color',[.2, .8, 1],'linew',2)
    fill(data_2_plot{6,4},squeeze(data_2_plot{6,3}(Combinations_2(cc,1),Combinations_2(cc,2),:)),[.2, .8, 1],'FaceAlpha',0.1,'EdgeColor','none')
    plot(freq_v,squeeze(data_2_plot{6,1}(Combinations_2(cc,1),Combinations_2(cc,2),:)),'Color',[.2, .8, 1],'linew',2)

    if          cc == 1
        legend('PL --> IL','', 'IL --> PL','')
        legend('boxoff')
        legend('FontSize',8)        
        ylabel({'ITI Freezing';'Granger causality'},'FontWeight','bold','FontSize',12)
    elseif      cc == 2
        legend('PL --> HPC','','HPC --> PL','')
        legend('boxoff')
        legend('FontSize',8)       
    else
        legend('IL --> HPC','', 'HPC --> IL','')
        legend('boxoff')
        legend('FontSize',8)      
    end

    xlabel('(Hz)','FontSize',12)
    xlim([2 11])
    ylim([0.005 6*10^-2])

end


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
