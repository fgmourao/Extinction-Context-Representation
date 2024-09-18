
figure
%plot(data.lfp{6,2}(3,1500:4500),'linew',1,'color', [0 0 0 .2])
hold on
plot(data.lfp{6,5}(1,11500:15500)./sum(abs(data.lfp{6,2}(1,11500:15500))),'linew',2)
%plot(data.lfp{6,5}(2,5000:9500))
plot(data.lfp{6,5}(3,11500:15500)./sum(abs(data.lfp{6,2}(3,11500:15500))),'linew',2)
ylim([-1.5*10^-3 1.5*10^-3])

figure
%plot(data.lfp{9,2}(3,2000:6500),'linew',1,'color', [0 0 0 .2])
hold on
plot(data.lfp{9,5}(1,:,5)/sum(abs(data.lfp{6,2}(1,6500:10500))),'linew',2);
%plot(data.lfp{9,5}(2,:,1),'linew',2);
plot(data.lfp{9,5}(3,:,5)/sum(abs(data.lfp{6,2}(1,6500:10500))),'linew',2);
ylim([-1*10^-3 1*10^-3])



figure


%plot(data.lfp{5,5}(3,data.events_behavior{2, 1}(2,1):data.events_behavior{2, 1}(2,2)),'linew',1)
plot(data.lfp{5,5}(3,:),'linew',1)


hold on
plot(data.lfp{9,5}(1,:,2),'linew',2);
%plot(data.lfp{9,5}(2,:,2),'linew',2);
plot(data.lfp{9,5}(3,:,2),'linew',2);
ylim([-300 300])





%% Plot to check all channels separately

figure
% Choose data set
data2plot_LFP = 5;

% Choose filter band
ff_1 = 5;
%ff_2 = 3;

% Choose behavior data set
data2plot_behavior = data.behavior_bins{1,1};

% Events
% CS indexes
cs_trial = round(data.events{2, 1}./parameters.downsampling);
% Freezing indexes
freezing_start = data.behavior{2, 1}(1,:);
freezing_end   = data.behavior{2, 1}(1,:)+data.behavior{2, 1}(2,:)-1;
% Time vector
behav_bins_time_v = linspace(1,length(data.behavior_bins{1,1})/(1/parameters.time_behav_bins),length(data.behavior_bins{1,1}));
behav_time_v = linspace(0,length(data.behavior{1,1})/1000,length(data.behavior{1,1}));

%yyaxis left
plot(data.timev_decimated,data.lfp{data2plot_LFP,ff_1}(1,:)./sum(abs(data.lfp{data2plot_LFP,ff_1}(1,:))))%,'Color','[0.3, 0.3, 0.3]')
hold on
plot(data.timev_decimated,data.lfp{data2plot_LFP,ff_1}(3,:)./sum(abs(data.lfp{data2plot_LFP,ff_1}(3,:))))%,'Color','[0.3, 0.3, 0.3]')

plot([data.timev_decimated(cs_trial(:,1));data.timev_decimated(cs_trial(:,2))], 10^-6.*([ones(1,length(cs_trial(:,1))).*1.5;ones(1,length(cs_trial(:,2))).*1.5]),'k-','linew', 20,'Color',[1, .4, .4])
%xline([data.timev_decimated(cs_trial(:,1))],'--r') % --> Exposure

ylabel('uV')

%plot(data.timev,data.lfp{data2plot_LFP,ff_2}(channels{jj}(ii),:),'Color','[0.6350, 0.0780, 0.1840]','linew',2)
%ylim([-600 600])

yyaxis right
%plot(behav_bins_time_v, movmean(data2plot_behavior,(1/parameters.time_behav_bins))-10^-5,'linew', 1, 'Color',[0.6, 0.6, 0.6, 0.6])

a = gca; % Get axis
a.YColor = 'w';
a.YTick = [];
hold all
plot([behav_time_v(freezing_start);behav_time_v(freezing_end)], [ones(1,length(freezing_start)).*-525;ones(1,length(freezing_end)).*-525],'k-','linew', 2,'Color',[.6, 0, 0])
%ylim([-500 0])
%         if ii == 1
%             legend('Signal - lowpass < 300 Hz','movement','Freezing','CS-Trials','location','northoutside','NumColumns',4)
%             legend('boxoff')
%         end

%ylim([min(data.data{1,ff}(ii,:)) max(data.data{1,ff}(ii,:))])
%ylim([-1000 1000])
%xlim([0 370]);
x1 = data.events_behavior{2, 1}(5,1)
x2 = data.events_behavior{2, 1}(5,2)

xlim([(x1/1000)-3 (x1/1000) + 3]);
