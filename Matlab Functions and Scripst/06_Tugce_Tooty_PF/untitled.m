%% Prepare all behavior data

data_fp_filt_lowpass = [];
data_fp = [];

for ii = 1:size(behavior,2)
    data_fp_filt_lowpass(ii,:) = eegfilt(Fiber_data{1, ii}(1:59400,2)',30,0,0.01);
    data_fp(ii,:) = Fiber_data{1, ii}(1:59400,2)';

end

data_fp_mean = mean(data_fp,1);
data_fp_filt_lowpass_mean = mean(data_fp_filt_lowpass,1);
data_fp_filt_lowpass_SEM  = std(data_fp_filt_lowpass,[],1)./size(data_fp_filt_lowpass,1);


data_fp_filt_highpass = [];

for ii = 1:size(behavior,2)
    data_fp_filt_highpass(ii,:) = eegfilt(Fiber_data{1, ii}(1:59400,2)',30,.01,[]);
end

data_fp_filt_highpass_mean = mean(data_fp_filt_highpass,1);

%% Prepare all behavior data

data_behavior = [];
data_behavior_filt = [];

for ii = 1:size(behavior,2)

    data_behavior(ii,:)      = behavior{1, ii}.data.behavior{1, 1}(1:59400);
    data_behavior_filt(ii,:) = eegfilt(behavior{1, ii}.data.behavior{1, 1}(1:59400),30,[],0.01);


end

data_behavior_mean = mean(data_behavior,1);
data_behavior_filt_mean = mean(data_behavior_filt,1);
data_behavior_filt_SEM = std(data_behavior_filt,[],1)./size(data_behavior_filt,1);


%% Select data to plot

figure
set(gcf,'color','w');
sc = [1,1,960,1200];
set(gcf, 'Position', sc);



subplot(4,4,[1 3 5 7 9 11 13 15])
% Movement
data_2_plot_mean_1       = (data_behavior_mean - min(data_behavior_mean))./max(data_behavior_mean);

% Delta F/F
data_2_plot_mean_2       = (data_fp_mean - min(data_fp_mean))./max(data_fp_mean);


% Time vector
behav_time_v = linspace(0,length(data_2_plot_mean_1)./behavior{1, 1}.parameters.original_srate,length(data_2_plot_mean_1));

% CS indexes
cs_trial       = behavior{1, 1}.data.events{1, 1};

% Freezing indexes
% freezing_ = cat(1,data.events_behavior{4,:});
% freezing_start_CS = freezing_(:,1);
% freezing_end_CS   = freezing_(:,2);
%
% freezing_ = cat(1,data.events_behavior{6,:});
% freezing_start_ITI = freezing_(:,1);
% freezing_end_ITI   = freezing_(:,2);

% Plot

%plot(behav_bins_time_v,data_2_plot,'Color',[0.3, 0.3, 0.3, 0.3]) % raw data
%plot(behav_time_v, data_2_plot_mean_1,'linew', 1, 'Color',[0.9, 0.9, 0.9, .5])

%plot(behav_time_v,ones(1,length(behav_time_v)).*200,'k--')
hold on
plot(behav_time_v, data_2_plot_mean_1_1+std(data_2_plot_mean_1,[],2),'linew', 3.5,'Color',[0.3, 0.3, 0.3, .5])

% plot([behav_time_v(freezing_start_CS);behav_time_v(freezing_end_CS)], [ones(1,length(freezing_start_CS)).*.85;ones(1,length(freezing_end_CS)).*.85],'k-','linew', 5,'Color',[.6, 0, 0,.8])
% plot([behav_time_v(freezing_start_ITI);behav_time_v(freezing_end_ITI)], [ones(1,length(freezing_start_ITI)).*.85;ones(1,length(freezing_end_ITI)).*.85],'k-','linew', 5,'Color',	[0, 0.4470, 0.7410, .5])

plot([behav_time_v(cs_trial(:,1));behav_time_v(cs_trial(:,2))], [ones(1,length(cs_trial(:,1))).*.98;ones(1,length(cs_trial(:,2))).*.98],'k-','linew', 20,'Color',[1, .4, .4])

xlabel('Time (s)','FontSize',12), ylabel('Movement ','FontSize',12)
xlim([behav_time_v(1)-10 behav_time_v(end)])
ylim([min(data_2_plot_mean_1) max(data_2_plot_mean_1)])


% Delta F/F
yyaxis right
%plot(behav_time_v,data_2_plot_mean_2,'color',[0 .6 1 .5],'linew', 1, 'LineStyle','-')
% hold on
plot(behav_time_v,data_2_plot_mean_2_2+std(data_2_plot_mean_2,[],2),'color',[0 .6 1],'linew', 3.5, 'LineStyle','-')

%yline(1.69,'k--')
%ylim([-2 6])
ylim([0 1])
% xlim([0 570])
ylabel('\Delta F/F ')



idx_to_plot = 4:4:28;

for ii = 1:length(idx_to_plot)

    subplot(7,4,idx_to_plot(ii))

    % Movement
    data_2_plot_1     = normalize(data_behavior(ii,:),'range');
    data_2_plot_1_1   = normalize(data_behavior_filt(ii,:),'range');

    % Delta F/F
    data_2_plot_2      = normalize(data_fp(ii,:),'range');
    data_2_plot_2_2    = normalize(data_fp_filt_lowpass(ii,:),'range')';

    %plot(behav_bins_time_v,data_2_plot,'Color',[0.3, 0.3, 0.3, 0.3]) % raw data
    % plot(behav_time_v, data_2_plot_1,'linew', 1, 'Color',[0.9, 0.9, 0.9])
    %plot(behav_time_v,ones(1,length(behav_time_v)).*200,'k--')
    hold on
    plot(behav_time_v, data_2_plot_1_1,'linew', 2,'Color',[0.3, 0.3, 0.3, .5])
    % plot([behav_time_v(freezing_start_CS);behav_time_v(freezing_end_CS)], [ones(1,length(freezing_start_CS)).*.85;ones(1,length(freezing_end_CS)).*.85],'k-','linew', 5,'Color',[.6, 0, 0,.8])
    % plot([behav_time_v(freezing_start_ITI);behav_time_v(freezing_end_ITI)], [ones(1,length(freezing_start_ITI)).*.85;ones(1,length(freezing_end_ITI)).*.85],'k-','linew', 5,'Color',	[0, 0.4470, 0.7410, .5])

    plot([behav_time_v(cs_trial(:,1));behav_time_v(cs_trial(:,2))], [ones(1,length(cs_trial(:,1))).*.95;ones(1,length(cs_trial(:,2))).*.95],'k-','linew', 20,'Color',[1, .4, .4])
    %ylim([min(data_2_plot_mean_1) 2500])
    ylim([0 1])

    xlabel('Time (s)','FontSize',12), ylabel('Movement ','FontSize',12)
    xlim([behav_time_v(1)-10 behav_time_v(end)])


    % Delta F/F
    yyaxis right
    %plot(behav_time_v,data_2_plot_2,'color',[0 .6 1 .3],'linew', 1, 'LineStyle','-')
%     hold on
    plot(behav_time_v,data_2_plot_2_2,'color',[0 .6 1],'linew', 2, 'LineStyle','-')

    %yline(1.69,'k--')
    ylim([0 1])
    % xlim([0 570])
    ylabel('\Delta F/F ')

end



