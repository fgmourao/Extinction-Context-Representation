% Cross correlation final plot

% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  01/2024
% Last update: 07/2024


%% Prepare data


%% Load data Extinction

% data_2_plot_x_corr_{1,ms} = x_corr.result_full;
% data_2_plot_x_corr_{2,ms} = x_corr.result;
% data_2_plot_x_corr_{3,ms} = x_corr.result_full_first_trials_mean;
% data_2_plot_x_corr_{4,ms} = x_corr.result_full_last_trials_mean;

for ii = 1:size(x_corr.result,1)
    for jj = 1:size(x_corr.result,2)

        result_full{ii,jj}(:,:,:,:,ms)                    = x_corr.result_full{ii,jj};
        result{ii,jj}(:,:,:,:,ms)                         = x_corr.result{ii,jj},;
        result_full_first_trials_mean{ii,jj}(:,:,:,:,ms)  = x_corr.result_full_first_trials_mean{ii,jj};
        result_full_last_trials_mean{ii,jj}(:,:,:,:,ms)   = x_corr.result_full_last_trials_mean{ii,jj};
        
        result_full_first_peak{ii,jj}(:,:,ms)                      = x_corr.result_full_first_trials_peak{ii,jj};
        result_full_last_peak{ii,jj}(:,:,ms)                       = x_corr.result_full_last_trials_peak{ii,jj};

        result_full_first_peak_overtime_all_animals{ii,jj}(:,:,ms)  = x_corr.result_full_first_trials_peak_overtime{ii,jj};
        result_full_last_peak_overtime_all_animals{ii,jj}(:,:,ms)   = x_corr.result_full_last_trials_peak_overtime{ii,jj};


    end
end

clear('ii','jj')



%%
if ms == length(files.FilesLoaded{1, 1})

    %% Frequency vector
    lags = x_corr.lags;

    % mean and SEM
for ii = 1:size(x_corr.result,1)
    for jj = 1:size(x_corr.result,2)

    data_2_plot_result_full_mean_animals                     = cellfun(@(x)mean(x,5),result_full,'UniformOutput',false);
    data_2_plot_result_full_SEM_animals                      = cellfun(@(x)std(x,[],5)./size(x,5),result_full,'UniformOutput',false);
    
    data_2_plot_result_mean_animals                          = cellfun(@(x)mean(x,5),result,'UniformOutput',false);
    data_2_plot_result_full_first_trials_mean_mean_animals   = cellfun(@(x)mean(x,5),result_full_first_trials_mean,'UniformOutput',false);
    data_2_plot_result_full_first_trials_mean_SEM_animals    = cellfun(@(x)std(x,[],5)./size(x,5),result_full_first_trials_mean,'UniformOutput',false);

    data_2_plot_result_full_last_trials_mean_mean_animals    = cellfun(@(x)mean(x,5),result_full_last_trials_mean,'UniformOutput',false);
    data_2_plot_result_full_last_trials_mean_SEM_animals     = cellfun(@(x)std(x,[],5)./size(x,5),result_full_last_trials_mean,'UniformOutput',false);

    end
end




for ii = 1:size(result_full,1)
    for jj = 1:size(result_full,2)
        
        if jj == 1 % baseline condition
            [~,data_2_plot_result_mean_animals_peak_idx{jj,ii}]     = max(data_2_plot_result_full_mean_animals{jj,ii},[],2);
            data_2_plot_result_mean_animals_peak{jj,ii}             = lags(data_2_plot_result_mean_animals_peak_idx{jj,ii});

        else
            [~,data_2_plot_result_mean_animals_first_trial_peak_idx{jj,ii}] = max(data_2_plot_result_full_first_trials_mean_mean_animals{jj,ii},[],2);
            data_2_plot_result_mean_animals_first_trial_peak{jj,ii}         = lags(data_2_plot_result_mean_animals_first_trial_peak_idx{jj,ii});
            
            [~,data_2_plot_result_mean_animals_last_trial_peak_idx{jj,ii}]  = max(data_2_plot_result_full_last_trials_mean_mean_animals{jj,ii},[],2);
            data_2_plot_result_mean_animals_last_trial_peak{jj,ii}          = lags(data_2_plot_result_mean_animals_last_trial_peak_idx{jj,ii});            
        end
    end
end



for ii = 1:size(result_full,1)
    for jj = 1:size(result_full,2)

        if jj == 1 % baseline condition
            data_2_plot_result_full_first_trials_overtime{jj,ii}              =  normalize(data_2_plot_result_mean_animals{jj,ii},3,'range');
            [~,data_2_plot_result_full_first_trials_peak_overtime_idx{jj,ii}] =  max(data_2_plot_result_full_first_trials_overtime{jj,ii},[],3);
            data_2_plot_result_full_first_trials_peak_overtime{jj,ii}         =  lags(data_2_plot_result_full_first_trials_peak_overtime_idx{jj,ii});

        else
            data_2_plot_result_full_first_trials_overtime{jj,ii}              =  normalize(mean(data_2_plot_result_mean_animals{jj,ii}(:,:,:,1:10),4),3,'range');
            [~,data_2_plot_result_full_first_trials_peak_overtime_idx{jj,ii}] =  max(data_2_plot_result_full_first_trials_overtime{jj,ii},[],3);
            data_2_plot_result_full_first_trials_peak_overtime{jj,ii}         =  lags(data_2_plot_result_full_first_trials_peak_overtime_idx{jj,ii});
            data_2_plot_result_full_first_trials_peak_overtime_mean{jj,ii}    =  mean(data_2_plot_result_full_first_trials_peak_overtime{jj,ii},2);

            data_2_plot_result_full_last_trials_overtime{jj,ii}               =  normalize(mean(data_2_plot_result_mean_animals{jj,ii}(:,:,:,11:20),4),3,'range');
            [~,data_2_plot_result_full_last_trials_peak_overtime_idx{jj,ii}]  =  max(data_2_plot_result_full_last_trials_overtime{jj,ii},[],3);
            data_2_plot_result_full_last_trials_peak_overtime{jj,ii}          =  lags(data_2_plot_result_full_last_trials_peak_overtime_idx{jj,ii});
            data_2_plot_result_full_last_trials_peak_overtime_mean{jj,ii}   =  mean(data_2_plot_result_full_last_trials_peak_overtime{jj,ii},2);

        end
    end
end

%% Plot the results

% choose filter
ff = 3;

%xlim([-200 200])

for jj = 2:3

    figure(jj+3);
    set(gcf,'color','w');
    sc = [1,1,800,1200];
    set(gcf, 'Position', sc);

    title_ = [{'mPFC PL lead <-----> mPFC IL lead'};{'mPFC PL lead <-----> dHPC lead'};{'mPFC IL lead <-----> dHPC lead'}];
    sgtitle([{'Cross-Correlation Over Time '};title_{jj}],'FontSize',14);


    % Baseline
    subplot(5,3,1);
    boundedline(lags,data_2_plot_result_full_mean_animals{1,ff}(jj,:),data_2_plot_result_full_SEM_animals{1,ff}(jj,:),'LineWidth',2,'color',[.6 .6 .6],'transparency',.4);
    hold on
    plot(lags(data_2_plot_result_mean_animals_peak_idx{1,ff}(jj,1)),data_2_plot_result_full_mean_animals{1,ff}(jj,data_2_plot_result_mean_animals_peak_idx{1,ff}(jj,1)),'ro','MarkerEdgeColor','r','MarkerFaceColor','r');
    title({['lag = ', num2str(data_2_plot_result_mean_animals_peak{1,ff}(1,jj)),' ms'];[]})
    ylabel([{'Normalized'};{'Crosscorrelation'}]);
    xlim([-200 200])
    %ylim([min(data_2_plot_result_full_mean_animals{1,ff}(jj,:)) max(data_2_plot_result_full_mean_animals{1,ff}(jj,:))])
    ylim([0 1])
    xline(0,'k--','LineWidth', 2)
    box off
    % legend('','peak lag','FontSize',10,'location', 'east')
    % legend('boxoff')


    subplot(5,3,4)
    histogram(data_2_plot_result_full_first_trials_peak_overtime{1,ff}(jj,:),'BinWidth',50,'FaceColor',[.6 .6 .6]);
    %histogram(squeeze(result_full_first_peak{1,ff}(jj,:)),'BinWidth',50,'FaceColor',[.6 .6 .6]);
    %histogram(squeeze(result_full_first_peak_overtime_all_animals{2,1}(jj,:,:)),'BinWidth',50,'FaceColor',[.6 .6 .6]);

    ylabel('Count');
    % a = fitdist(xcorr_result_lag_peak','normal');
    % hold on
    % plot(a)
    xlim([-200 200])
    box off
    %ylim([0 10])
    xline(0,'k--','LineWidth', 2)


    subplot(5,3,[7 10 13])
    imagesc(lags, (1:size(data_2_plot_result_full_first_trials_overtime{1,ff},2)), normalize(squeeze(data_2_plot_result_full_first_trials_overtime{1,ff}(jj,:,:)),2,'range'));
    %contourf(lags,(1:size(x_corr.result{1,ff},2)),squeeze(x_corr.result{1,ff}(jj,:,:)),80,'linecolor','none');
    xlim([-200 200])
    set(gca, 'YDir','normal')
    xline(0,'w--','LineWidth', 2)
    c = colorbar('Ticks',[0.01,0.99],'TickLabels',{'min','max'},'location','southoutside');
    %c.Label.String = 'Normalized peak';
    colormap('parula')
    % Color map matplotlib
    % py_path = "~/anaconda3/envs/Python_3_10/bin/python";
    %     Py_map = getPyPlot_cMap('parula', [], [], py_path);
    %     colormap(Py_map)
    xlabel('Lag (ms)');
    ylabel('Time (seconds)');

    yyaxis right
    plot(data_2_plot_result_full_first_trials_peak_overtime{1,ff}(jj,:),(1:size(data_2_plot_result_full_first_trials_overtime{1,ff},2)),'wo','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',4);
    xlim([-200 200])
    ylim([0.5 size(x_corr.result{1,ff},2)+.5])
    a = gca; % Get axis
    a.YColor = 'w';
    a.YTick = [];
    % legend(r,'peak','location','southoutside','FontSize',14)
    % legend('boxoff')

    % CS-Tone. First 10 trials

    subplot(5,3,2)
    boundedline(lags,data_2_plot_result_full_first_trials_mean_mean_animals{2,ff}(jj,:),data_2_plot_result_full_first_trials_mean_SEM_animals{2,ff}(jj,:),'LineWidth',2,'color',[.6 .6 .6],'transparency',.4);
    hold on
    plot(lags(data_2_plot_result_mean_animals_first_trial_peak_idx{2,ff}(jj,1)),data_2_plot_result_full_first_trials_mean_mean_animals{2,ff}(jj,data_2_plot_result_mean_animals_first_trial_peak_idx{2,ff}(jj,1)),'ro','MarkerEdgeColor','r','MarkerFaceColor','r');
    title({['lag = ', num2str(data_2_plot_result_mean_animals_first_trial_peak{2,ff}(1,jj)),' ms'];[]})
    ylabel([{'Normalized'};{'Crosscorrelation'}]);
    xlim([-200 200])
    %ylim([min(data_2_plot_result_full_mean_animals_first_trials_mean {2,ff}(jj,:)) max(data_2_plot_result_full_mean_animals_first_trials_mean {2,ff}(jj,:))])
     ylim([0 1])
    xline(0,'k--','LineWidth', 2)
    box off
    % legend('','peak lag','FontSize',10,'location', 'east')
    % legend('boxoff')

    subplot(5,3,5)
    histogram(data_2_plot_result_full_first_trials_peak_overtime{2,ff}(jj,:),'BinWidth',50,'FaceColor',[.6 .6 .6]);
    %histogram(squeeze(result_full_first_peak{2,ff}(jj,:,:)),'BinWidth',50,'FaceColor',[.6 .6 .6]);
    %histogram(result_full_first_peak_overtime_all_animals{2,1}(jj,:,:),'BinWidth',50,'FaceColor',[.6 .6 .6]);

    ylabel('Count');
    % a = fitdist(xcorr_result_lag_peak','normal');
    % hold on
    % plot(a)
    xlim([-200 200])
    box off
    %ylim([0 10])
    xline(0,'k--','LineWidth', 2)

    subplot(5,3,[8 11 14])
    imagesc(lags, (1:size(data_2_plot_result_full_first_trials_overtime{2,ff},2)), normalize(squeeze(data_2_plot_result_full_first_trials_overtime{2,ff}(jj,:,:)),2,'range'));

    %contourf(lags,(1:size(x_corr.result{2,ff},2)), mean(squeeze(x_corr.result{2,ff}(jj,:,:,1:size(CSIT{ms},2)/2)),3),80,'linecolor','none');

    xlim([-200 200])
    set(gca, 'YDir','normal')
    xline(0,'w--','LineWidth', 2)
    c = colorbar('Ticks',[0.01,0.99],'TickLabels',{'min','max'},'location','southoutside');
    %c.Label.String = 'Normalized peak';
    colormap('parula')
    % Color map matplotlib
    % py_path = "~/anaconda3/envs/Python_3_10/bin/python";
    %     Py_map = getPyPlot_cMap('parula', [], [], py_path);
    %     colormap(Py_map)
    xlabel('Lag (ms)');
    ylabel('Time (seconds)');

    yyaxis right
    plot(data_2_plot_result_full_first_trials_peak_overtime{2,ff}(jj,:),(1:size(data_2_plot_result_full_first_trials_overtime{2,ff},2)),'wo','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',4);
    xlim([-200 200])
    ylim([0.5 size(x_corr.result{2,ff},2)+.5])
    a = gca; % Get axis
    a.YColor = 'w';
    a.YTick = [];
    % legend(r,'peak','location','southoutside','FontSize',14)
    % legend('boxoff')


% CS-Tone. last 10 trials

    subplot(5,3,3)
    boundedline(lags,data_2_plot_result_full_last_trials_mean_mean_animals{2,ff}(jj,:),data_2_plot_result_full_last_trials_mean_SEM_animals{2,ff}(jj,:),'LineWidth',2,'color',[.6 .6 .6],'transparency',.4);
    hold on
    plot(lags(data_2_plot_result_mean_animals_last_trial_peak_idx{2,ff}(jj,1)),data_2_plot_result_full_last_trials_mean_mean_animals{2,ff}(jj,data_2_plot_result_mean_animals_last_trial_peak_idx{2,ff}(jj,1)),'ro','MarkerEdgeColor','r','MarkerFaceColor','r');
    title({['lag = ', num2str(data_2_plot_result_mean_animals_last_trial_peak{2,ff}(1,jj)),' ms'];[]})
    ylabel([{'Normalized'};{'Crosscorrelation'}]);
    xlim([-200 200])
    %ylim([min(data_2_plot_result_full_mean_animals_last_trials_mean {2,ff}(jj,:)) max(data_2_plot_result_full_mean_animals_last_trials_mean {2,ff}(jj,:))])
     ylim([0 1])
    xline(0,'k--','LineWidth', 2)
    box off
    % legend('','peak lag','FontSize',10,'location', 'east')
    % legend('boxoff')

    subplot(5,3,6)
    histogram(data_2_plot_result_full_last_trials_peak_overtime{2,ff}(jj,:),'BinWidth',50,'FaceColor',[.6 .6 .6]);
    %histogram(squeeze(result_full_last_peak{2,ff}(jj,:,:)),'BinWidth',50,'FaceColor',[.6 .6 .6]);
    %histogram(result_full_last_peak_overtime_all_animals{2,1}(jj,:,:),'BinWidth',50,'FaceColor',[.6 .6 .6]);

    ylabel('Count');
    % a = fitdist(xcorr_result_lag_peak','normal');
    % hold on
    % plot(a)
    xlim([-200 200])
    box off
    %ylim([0 10])
    xline(0,'k--','LineWidth', 2)

    subplot(5,3,[9 12 15])
    imagesc(lags, (1:size(data_2_plot_result_full_last_trials_overtime{2,ff},2)), normalize(squeeze(data_2_plot_result_full_last_trials_overtime{2,ff}(jj,:,:)),2,'range'));

    %contourf(lags,(1:size(x_corr.result{2,ff},2)), mean(squeeze(x_corr.result{2,ff}(jj,:,:,1:size(CSIT{ms},2)/2)),3),80,'linecolor','none');

    xlim([-200 200])
    set(gca, 'YDir','normal')
    xline(0,'w--','LineWidth', 2)
    c = colorbar('Ticks',[0.01,0.99],'TickLabels',{'min','max'},'location','southoutside');
    %c.Label.String = 'Normalized peak';
    colormap('parula')
    % Color map matplotlib
    % py_path = "~/anaconda3/envs/Python_3_10/bin/python";
    %     Py_map = getPyPlot_cMap('parula', [], [], py_path);
    %     colormap(Py_map)
    xlabel('Lag (ms)');
    ylabel('Time (seconds)');

    yyaxis right
    plot(data_2_plot_result_full_last_trials_peak_overtime{2,ff}(jj,:),(1:size(data_2_plot_result_full_last_trials_overtime{2,ff},2)),'wo','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',4);
    xlim([-200 200])
    ylim([0.5 size(x_corr.result{2,ff},2)+.5])
    a = gca; % Get axis
    a.YColor = 'w';
    a.YTick = [];
    % legend(r,'peak','location','southoutside','FontSize',14)
    % legend('boxoff')

end

end

%% Save Figures

% path = '/Users/flavio/Desktop';
% %path = files.FilesLoaded{1, 1}.folder;
% 
% name = strcat(path,'/','_x_corr_extinction_6_8Hz_mean');
% 
% %get all figures
% figlist = findobj(allchild(0), 'flat', 'Type', 'figure');
% 
% %name_figs = {'_3_5Hz','_6_8Hz','_7_10Hz'};
% name_figs = {'_PL_IL','_PL_dHPC','_IL_dHPC'};
% 
% set(gcf,'renderer','Painters')
% 
% % Loop through figure
% for ii = 1:numel(figlist)
%     name_loop = strcat(name,name_figs{ii});
%     FigHandle = figlist(ii);
%     saveas(FigHandle,name_loop,'png')
%     exportgraphics(gcf,strcat(name_loop,'.eps'),'Resolution', 300)
% 
% end
% 
% close all
% clear('name','newStr','path')