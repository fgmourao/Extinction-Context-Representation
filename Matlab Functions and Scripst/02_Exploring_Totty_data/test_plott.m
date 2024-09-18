%% Plot the results

% choose filter
ff = 2;

%xlim([-x_corr.parameters.nlags x_corr.parameters.nlags])

for jj = 1:3

    figure(jj);
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
    xlim([-x_corr.parameters.nlags x_corr.parameters.nlags])
    %ylim([min(data_2_plot_result_full_mean_animals{1,ff}(jj,:)) max(data_2_plot_result_full_mean_animals{1,ff}(jj,:))])
    ylim([0.2 .5])
    xline(0,'k--','LineWidth', 2)
    box off
    % legend('','peak lag','FontSize',10,'location', 'east')
    % legend('boxoff')


    subplot(5,3,4)
    histogram(data_2_plot_result_full_first_trials_peak_overtime{1,ff}(jj,:),'BinWidth',20,'FaceColor',[.6 .6 .6]);
    ylabel('Count');
    % a = fitdist(xcorr_result_lag_peak','normal');
    % hold on
    % plot(a)
    xlim([-x_corr.parameters.nlags x_corr.parameters.nlags])
    box off
    %ylim([0 10])
    xline(0,'k--','LineWidth', 2)


    subplot(5,3,[7 10 13])
    imagesc(lags, (1:size(data_2_plot_result_full_first_trials_overtime{1,ff},2)), normalize(squeeze(data_2_plot_result_full_first_trials_overtime{1,ff}(jj,:,:)),2,'range'));
    %contourf(lags,(1:size(x_corr.result{1,ff},2)),squeeze(x_corr.result{1,ff}(jj,:,:)),80,'linecolor','none');
    xlim([-x_corr.parameters.nlags x_corr.parameters.nlags])
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
    xlim([-x_corr.parameters.nlags x_corr.parameters.nlags])
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
    xlim([-x_corr.parameters.nlags x_corr.parameters.nlags])
    %ylim([min(data_2_plot_result_full_mean_animals_first_trials_mean {2,ff}(jj,:)) max(data_2_plot_result_full_mean_animals_first_trials_mean {2,ff}(jj,:))])
    ylim([.4 .65])
    xline(0,'k--','LineWidth', 2)
    box off
    % legend('','peak lag','FontSize',10,'location', 'east')
    % legend('boxoff')

    subplot(5,3,5)
    histogram(data_2_plot_result_full_first_trials_peak_overtime{2,ff}(jj,:),'BinWidth',30,'FaceColor',[.6 .6 .6]);
    ylabel('Count');
    % a = fitdist(xcorr_result_lag_peak','normal');
    % hold on
    % plot(a)
    xlim([-x_corr.parameters.nlags x_corr.parameters.nlags])
    box off
    %ylim([0 10])
    xline(0,'k--','LineWidth', 2)

    subplot(5,3,[8 11 14])
    imagesc(lags, (1:size(data_2_plot_result_full_first_trials_overtime{2,ff},2)), normalize(squeeze(data_2_plot_result_full_first_trials_overtime{2,ff}(jj,:,:)),2,'range'));

    %contourf(lags,(1:size(x_corr.result{2,ff},2)), mean(squeeze(x_corr.result{2,ff}(jj,:,:,1:size(CSIT{ms},2)/2)),3),80,'linecolor','none');

    xlim([-x_corr.parameters.nlags x_corr.parameters.nlags])
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
    xlim([-x_corr.parameters.nlags x_corr.parameters.nlags])
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
    xlim([-x_corr.parameters.nlags x_corr.parameters.nlags])
    %ylim([min(data_2_plot_result_full_mean_animals_last_trials_mean {2,ff}(jj,:)) max(data_2_plot_result_full_mean_animals_last_trials_mean {2,ff}(jj,:))])
    ylim([.4 .65])
    xline(0,'k--','LineWidth', 2)
    box off
    % legend('','peak lag','FontSize',10,'location', 'east')
    % legend('boxoff')

    subplot(5,3,6)
    histogram(data_2_plot_result_full_last_trials_peak_overtime{2,ff}(jj,:),'BinWidth',30,'FaceColor',[.6 .6 .6]);
    ylabel('Count');
    % a = fitdist(xcorr_result_lag_peak','normal');
    % hold on
    % plot(a)
    xlim([-x_corr.parameters.nlags x_corr.parameters.nlags])
    box off
    %ylim([0 10])
    xline(0,'k--','LineWidth', 2)

    subplot(5,3,[9 12 15])
    imagesc(lags, (1:size(data_2_plot_result_full_last_trials_overtime{2,ff},2)), normalize(squeeze(data_2_plot_result_full_last_trials_overtime{2,ff}(jj,:,:)),2,'range'));

    %contourf(lags,(1:size(x_corr.result{2,ff},2)), mean(squeeze(x_corr.result{2,ff}(jj,:,:,1:size(CSIT{ms},2)/2)),3),80,'linecolor','none');

    xlim([-x_corr.parameters.nlags x_corr.parameters.nlags])
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
    xlim([-x_corr.parameters.nlags x_corr.parameters.nlags])
    ylim([0.5 size(x_corr.result{2,ff},2)+.5])
    a = gca; % Get axis
    a.YColor = 'w';
    a.YTick = [];
    % legend(r,'peak','location','southoutside','FontSize',14)
    % legend('boxoff')

end
