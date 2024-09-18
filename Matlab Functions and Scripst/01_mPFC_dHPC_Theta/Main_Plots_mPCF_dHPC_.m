
%% Plots
%  Call functions/scripts

% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  12/2023
% Last update:

%% Check frequency band based on data cursor/cursor info

% f = 1/(cursor_info(1,1).Position(1,1)-cursor_info(1,2).Position(1,1));

%% Plot to check all channels separately

% Choose data set
data2plot_LFP = 5;

% Choose filter band
ff_1 = 1;
%ff_2 = 3;

% Choose behavior data set
data2plot_behavior = data.behavior_bins{1,1};

% Events
% CS indexes
cs_trial = round(data.events{2, 1}./parameters.downsampling); 
% Freezing indexes
freezing_start = data.behavior_bins{2, 1}(1,:);
freezing_end   = data.behavior_bins{2, 1}(1,:)+data.behavior_bins{2, 1}(2,:)-1;
% Time vector
behav_bins_time_v = linspace(0,length(data.behavior_bins{1,1})/(1/parameters.time_behav_bins),length(data.behavior_bins{1,1}));


% Choose channels to plot
% Totty data -> 1-7   mPFC PL
%            -> 8-14  mPFC IL
%            -> 15-29 HPC

channels{1} = 1;
channels{2} = 2;
channels{3} = 3;

figure
set(gcf,'color','w');
set(gcf, 'Position', get(0, 'Screensize'));

for jj = 1:length(channels)
    for ii = 1:length(channels{jj})

%         if jj == 1 || jj == 2 % Conditioning for several channels per substract
%             figure(jj)
%             set(gcf,'color','w');
%             set(gcf, 'Position', get(0, 'Screensize'));
%             subplot(3,4,ii)
% 
%         else
%             figure(jj)
%             set(gcf,'color','w');
%             set(gcf, 'Position', get(0, 'Screensize'));
%             subplot(4,4,ii)
%         end


        subplot(3,1,jj)

        %yyaxis left
        plot(data.timev_decimated,data.lfp{data2plot_LFP,ff_1}(channels{jj}(ii),:),'Color','[0.3, 0.3, 0.3]')
        hold all
        plot([data.timev_decimated(cs_trial(:,1));data.timev_decimated(cs_trial(:,2))], 500.*([ones(1,length(cs_trial(:,1))).*.5;ones(1,length(cs_trial(:,2))).*.5]),'k-','linew', 20,'Color',[1, .4, .4])  
        %xline([data.timev_decimated(cs_trial(:,1))],'--r') % --> Exposure 

        ylabel('uV')

        %plot(data.timev,data.lfp{data2plot_LFP,ff_2}(channels{jj}(ii),:),'Color','[0.6350, 0.0780, 0.1840]','linew',2)
        ylim([-500 500])

        yyaxis right
        plot(behav_bins_time_v, movmean(data2plot_behavior,(1/parameters.time_behav_bins))-500,'linew', 1, 'Color',[0.6, 0.6, 0.6, 0.6])

        a = gca; % Get axis
        a.YColor = 'w';
        a.YTick = [];
        hold all
        plot([behav_bins_time_v(freezing_start);behav_bins_time_v(freezing_end)], [ones(1,length(freezing_start)).*-425;ones(1,length(freezing_end)).*-425],'k-','linew', 2,'Color',[.6, 0, 0])
        ylim([-500 0])
%         if ii == 1
%             legend('Signal - lowpass < 300 Hz','movement','Freezing','CS-Trials','location','northoutside','NumColumns',4)
%             legend('boxoff')
%         end

        %ylim([min(data.data{1,ff}(ii,:)) max(data.data{1,ff}(ii,:))])
        %ylim([-1000 1000])
        %xlim([0 370]);
        xlim([0 data.timev_decimated(1,end)]);

        %xlim([0 20]);

        %yticklabels({'-1','0','1'})

        if jj == 1
            title(['mPFC IL ',num2str(ii)])
            a.TitleHorizontalAlignment = 'left';
        elseif jj == 2
            title(['mPFC PL ',num2str(ii)])
            a.TitleHorizontalAlignment = 'left';
        else
            title(['HPC ',num2str(ii)])
            a.TitleHorizontalAlignment = 'left';
        end

        xlabel('Time (s)')

        box off

        sgtitle({[id;'']})
    end

end

%% Save

newStr = id(1:end-5);

path = '/Users/flavio/Desktop';
%path = files.FilesLoaded{1, 1}.folder;

name = strcat(path,'/',newStr,'_main_plots');

saveas(gcf,name,'png')

close all

clear('name','newStr1','path')
clear ('data2plot_LFP','ff_1','data2plot_behavior','cs_trial','freezing_start','freezing_end','channels','jj','a','ii','newStr1','path')

%% Plot to check all channels in the same plot

% % Choose data set
% data2plot_LFP = 5;
% 
% % Choose filter band
% ff_1 = 2;
% %ff_2 = 2;
% 
% % Choose behavior data set
% data2plot_behavior = data.behavior_bins{1,1};
% 
% % Events
% % CS indexes
% cs_trial = round(data.events{2, 1});
% % Freezing indexes
% freezing_start = data.behavior_bins{2, 1}(1,:);
% freezing_end   = data.behavior_bins{2, 1}(1,:)+data.behavior_bins{2, 1}(2,:);
% % Time vector
% behav_bins_time_v = linspace(1,length(data.behavior_bins{1,1})/(1/parameters.time_behav_bins),length(data.behavior_bins{1,1}));
% 
% 
% % Choose channels to plot
% channels = 2:16;
% 
% % factor
% factor = (channels)'*1500;
% 
% % Set Figure
% figure
% set(gcf,'color','white')
% box 'off'
% hold on
% 
% % Select fields with data
% r = plot(data.timev_decimated, bsxfun(@plus, data.lfp{data2plot_LFP, ff_1}(channels,:), factor),'Color','[0.3, 0.3, 0.3]','linew',1);
% 
% a = gca;
% %a.YLim = [-500  3000];
% a.YColor = 'w';
% a.YTick = [];
% a.XLim = [0  data.timev_decimated(end)];
% 
% title('Habituation')
% a.TitleHorizontalAlignment = 'left';
% 
% yyaxis right
% plot(behav_bins_time_v, movmean(data2plot_behavior,(1/parameters.time_behav_bins))-400,'linew', 1, 'Color',[0.6, 0.6, 0.6, 0.6])
% 
% hold all
% plot([behav_bins_time_v(freezing_start);behav_bins_time_v(freezing_end)], [ones(1,length(freezing_start)).*-50;ones(1,length(freezing_end)).*-50],'k-','linew', 2,'Color',[.6, 0, 0])
% %plot([data.timev_decimated(cs_trial(:,1));data.timev_decimated(cs_trial(:,2))], [ones(1,length(cs_trial(:,1))).*-15;ones(1,length(cs_trial(:,2))).*-15],'k-','linew', 20,'Color',[1, .4, .4,.3])
% 
% % Plot sound epochs
% %I = plot([data.events.idx_t(:)';data.events.idx_t(:)'], [zeros(1,length(data.events.idx_t(:)));a.YLim(2)*(ones(1,length(data.events.idx_t(:))))],'Color',[0.6350, 0.0780, 0.1840]','linew',2);
% 
% a = gca; % Get axis
% a.YColor = 'w';
% a.YTick = [];
% box off
% 
% %a.XLim = [67.6  69.6];
% %a.XLim = [127.6  129.6];
% 
% % a.XLim = [ 91  93];
% % a.XLim = [181 183];
% % a.XLim = [301 303];
% % a.XLim = [391 393];
% % a.XLim = [511.4 513];
% 
% % a.XLim = [175.2 177.2];
% % a.XLim = [295.2 297.2];
% % a.XLim = [385.2 387.2];
% % a.XLim = [505.2 507.2];
% 
% % a.XLim = [295.2 297];
% % a.XLim = [235.2 237];
% % a.XLim = [175.2 177];
% % a.XLim = [115.2 117];
% % a.XLim = [55.2 77];
% 
% % a.XLim = [55 57];
% % a.XLim = [115 117];
% % a.XLim = [175.1 177.1];
% 
% % a.XLim = [295 297];
% 
% xlabel('Time (s)')
% 
% legend('Signal - lowpass < 300 Hz','movement','Freezing','CS-Trials','location','northoutside','NumColumns',4)
% legend('boxoff')
% 
% % Clear trash
% clear ('factor','channels','filter','substrate','session','str','sub','r','a','I','lh','lh_pos');

%% Pre Processing - Plot to check Trials
% 

% ----- CHANGE THIS CODE -----

% % Choose filter band
% ff1 = 4;
% ff2 = 3;
% 
% % Choose channel to plot
% ch = 13;
% 
% % Set Figure
% figure
% set(gcf,'color','white')
% 
% 
% for ii = 1:parameters.NTrials
%     subplot(1,5,ii)
%     hold on
%     plot(data.time_trials,data.data_trials{1,ff1}(ch,:,ii),'Color', '[0.7 0.7 0.7]','linew',2)
%     plot(data.time_trials,data.data_trials{1,ff2}(ch,:,ii),'Color', '[0.2, 0.6, 1.0]','linew',2)
% 
%     plot(data.time_trials,data.data_trials{1,ff2}(1,:,ii).* 20 + 300,'Color','[0.6350, 0.0780, 0.1840, 0.2]','linew',2)
%     
%     ylim([-500 500]);
%     xlim([-parameters.Tpre parameters.trialperiod + parameters.Tpos]);
%     title(['Trial ',num2str(ii)])
%     xlabel('Time (s)')
%     ylabel('mV')
%     box off
% 
% end
%   
% clear ('ii','ff','ch')

%% last update: 18/01/24
%  listening: 
