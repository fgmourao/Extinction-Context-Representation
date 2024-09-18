% Welch`s final plot

% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  01/2024
% Last update: 06/2024


%% Prepare data


%% Load data Extinction

%data_2_plot_extinction_mean_all = [];

for ii = 1:size(pw.full_trial.Pxx_TotalPower_norm,1)

    data_2_plot_extinction_mean_all{ii,ms} = pw.full_trial.Pxx_TotalPower_norm{ii,1};

end


%data_2_plot_extinction_mean_first_trials = [];

for ii = 1:size(pw.full_trial.Pxx_TotalPower_norm_mean,1)

    data_2_plot_extinction_mean_first_trials{ii,ms} = pw.full_trial.Pxx_TotalPower_norm_mean{ii,1};

end




if ms == length(files.FilesLoaded{1, 1})

    %% Frequency vector
    freq_v = pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz);

    %% Concatenate animals
    data_2_plot_extinction_mean_SEM{1,1} = cat(3,data_2_plot_extinction_mean_all{1,:});


    %% Mean and SEM

    data_2_plot_extinction_mean_all_SEM          = [];
    data_2_plot_extinction_mean_first_trials_SEM = [];
    data_2_plot_extinction_mean_last_trials_SEM  = [];


    for ii = 1:size(data_2_plot_extinction_mean_all,1)
        data_2_plot_extinction_mean_all_SEM{ii,1}  = mean(cat(4,data_2_plot_extinction_mean_all{ii,:}),4,'omitnan');
        data_2_plot_extinction_mean_all_SEM{ii,2}  = std(cat(4,data_2_plot_extinction_mean_all{ii,:}),[],4,'omitnan')./size(cat(4,data_2_plot_extinction_mean_all{ii,:}),4);

        data_2_plot_extinction_mean_first_trials_SEM{ii,1}  = mean(cat(3,data_2_plot_extinction_mean_first_trials{ii,:}),3,'omitnan');
        data_2_plot_extinction_mean_first_trials_SEM{ii,2}  = std(cat(3,data_2_plot_extinction_mean_first_trials{ii,:}),[],3,'omitnan')./size(cat(3,data_2_plot_extinction_mean_first_trials{ii,:}),3);

    end

    clear('ii')

    %% Plot to check - for now I am considring just CS-Tone

    %data_2_plot = pw.full_trial.Pxx_stats_baseline;
    data_2_plot       = data_2_plot_extinction_mean_all_SEM;
    data_2_plot_mean1 = data_2_plot_extinction_mean_first_trials_SEM;

    % choose channel
    ch1 = 1;
    ch2 = 2;
    ch3 = 3;


    f = figure;
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gcf,'color','w');
    axis('square')
    sgtitle('\fontsize{18} \bf Extinction')

    % Baseline
    b(1) = subplot(3,7,1);
    boundedline(freq_v,data_2_plot{1, 1}(ch1,:),data_2_plot{1, 2}(ch1,:),'linewidth',2,'color',[.6 .6 .6],'transparency',.4)
    xlim([2 12])
    ylabel({'\fontsize{16} \bf mPFC (PL)';'Normalized Amplitude (a.u.)';[]})
    xlabel('Hz')
    yticks([0 0.05 0.1 0.15 0.2])
    box off
    title({'\fontsize{16} \bf Baseline';[]})

    b(2) = subplot(3,7,8);
    boundedline(freq_v,data_2_plot{1, 1}(ch2,:),data_2_plot{1, 2}(ch2,:),'linewidth',2,'color',[.6 .6 .6],'transparency',.4 )
    xlim([2 12])
    ylabel({'\fontsize{16} \bf mPFC (IL)';'Normalized Amplitude (a.u.)';[]})
    box off
    xlabel('Hz')
    yticks([0 0.05 0.1 0.15 0.2])

    b(3) = subplot(3,7,15);
    boundedline(freq_v,data_2_plot{1, 1}(ch3,:),data_2_plot{1, 2}(ch3,:),'linewidth',2,'color',[.6 .6 .6],'transparency',.4 )
    xlim([2 12])
    box off
    ylabel({'\fontsize{16} \bf dHPC';'Normalized Amplitude (a.u.)';[]})
    xlabel('Hz')
    yticks([0 0.05 0.1 0.15 0.2])



    % CS-Tones

    for ii = 1:5

        % First 10 trials
        t(ii+1) = subplot(3,7,ii+1);
        boundedline(freq_v,data_2_plot{2, 1}(ch1,:,ii),data_2_plot{2, 2}(ch1,:,ii),'linewidth',2,'Color',[0 0 .6],'transparency',.4)

        xlim([2 12])
        %ylabel('Amplitude')
        xlabel('Hz')
        axis('square')
        box off
        %title({['CS-Tone ' num2str(CSIT{ms}(1,ii))];[]})

        t(ii+25) = subplot(3,7,ii+8);
        boundedline(freq_v,data_2_plot{2, 1}(ch2,:,ii),data_2_plot{2, 2}(ch2,:,ii),'linewidth',2,'Color',[0 0 .6],'transparency',.4)

        xlim([2 12])
        %ylabel('Amplitude')
        axis('square')
        box off
        xlabel('Hz')
        title({['CS-Tone ' num2str(ii)];[]})

        t(ii+49) = subplot(3,7,ii+15);
        boundedline(freq_v,data_2_plot{2, 1}(ch3,:,ii),data_2_plot{2, 2}(ch3,:,ii),'linewidth',2,'Color',[0 0 .6],'transparency',.4)

        xlim([2 12])
        %ylabel('Amplitude')
        axis('square')
        box off
        xlabel('Hz')
        title({['CS-Tone ' num2str(ii)];[]})

    end


    % Trial averaged

    m(1) = subplot(3,7,7);
    boundedline(freq_v,data_2_plot_mean1{2, 1}(ch1,:),data_2_plot_mean1{2, 2}(ch1,:),'linewidth',2,'color',[0 0 .6],'transparency',.4 )
 
    xlim([2 12])
    xlabel('Hz')
    yticks([0 0.05 0.1 0.15 0.2])
    box off
    title({'\fontsize{16} \bf Averaged Trials';[]})


    m(2) = subplot(3,7,14);
    boundedline(freq_v,data_2_plot_mean1{2, 1}(ch2,:),data_2_plot_mean1{2, 2}(ch2,:),'linewidth',2,'color',[0 0 .6],'transparency',.4 )
 
    xlim([2 12])
    box off
    xlabel('Hz')
    yticks([0 0.05 0.1 0.15 0.2])
    %title({'Averaged Trials'})


    m(3) = subplot(3,7,21);
    boundedline(freq_v,data_2_plot_mean1{2, 1}(ch3,:),data_2_plot_mean1{2, 2}(ch3,:),'linewidth',2,'color',[0 0 .6],'transparency',.4 )

    xlim([2 12])
    box off
    xlabel('Hz')
    yticks([0 0.05 0.1 0.15 0.2])
    %title({'Averaged Trials'})
    legend('\fontsize{14} First 10 Trials','\fontsize{14} Last 10 Trials')
    lg = legend;
    set(lg,...
        'Position',[0.8859375 0.290602184305253 0.0497395833333333 0.0241788321167883]);

    linkaxes([b t m],'xy');
    %b(1).YLim  = [0 15];
    b(1).YLim  = [0 0.10];


    clear ('ii','ch1','ch2','ch3','f','lg','m','t','b','data_2_plot','data_2_plot_mean1','data_2_plot_mean2')

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
