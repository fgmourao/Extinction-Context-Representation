
%% Short-time FFT by matlab built function spectrogram

% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  12/2023
% Last update: 05/2024

%%
% First extract the data with: Extracting_LFPs_and_events.m
% ... then organize the data with the script: Pre_processing.m

%%

fprintf('\n sFFT Spectrogram... \n');

%%

% Initializing
short_fft = [];

% Time window
short_fft.timewin    = 25000; % in ms

% Convert time window to points
short_fft.timewinpnts  = hamming(round(short_fft.timewin/(1000/parameters.decimated_srate)));

% nFFT
%short_fft.nFFT = 2^nextpow2(round(short_fft.timewin/(1000/parameters.decimated_srate)));
short_fft.nFFT = 2^16; 

% Number of overlap samples
short_fft.overlap = 90;
short_fft.noverlap = floor(short_fft.overlap*0.01*(round(short_fft.timewin/(1000/parameters.decimated_srate))));


% Spectrogram
% lines: frequencies / columns: time / third dimension: channels

% Choose data according pre_processing.m define
not = 5;

for ii = 1:size(data.lfp{not,1},1)
    if ii == 1
       [short_fft.data(:,:,ii),short_fft.freq,short_fft.time] = spectrogram(data.lfp{not,1}(ii,:),short_fft.timewinpnts,short_fft.noverlap,short_fft.nFFT,parameters.decimated_srate);
    else
        short_fft.data(:,:,ii) = spectrogram(data.lfp{not,1}(ii,:),short_fft.timewinpnts,short_fft.noverlap,short_fft.nFFT,parameters.decimated_srate);
    end
end

clear ('ii','jj','not')

%% Define indexes from spectrogram - Considering the entire trial period

% rows -> trials
% columns 1 -> CS_Trial star 
% columns 2 -> sound start
% columns 3 -> sound stop
% columns 4 -> time after sound epoch

% CS-Trials epochs
short_fft.time_idx_t = (data.events{3, 1}(:)./parameters.decimated_srate);
short_fft.time_idx = reshape(dsearchn(short_fft.time',short_fft.time_idx_t),5,2);
short_fft.time_idx_t = reshape(short_fft.time_idx_t,5,2);  % just to keep the same format m x n


%% Plot to check full session. Channels per substrate 

% Choose channel
ch = 1:3;

%Define frequencies to plot in each subplot
steps = diff(short_fft.freq); % according to the fft time window

% For data
freq2plot = 3:steps(1):9;
closestfreq = dsearchn(short_fft.freq,freq2plot');

% For events
freq2plot_events = 1:steps(1):5;
closestfreq_events = dsearchn(short_fft.freq,freq2plot_events');

% Choose behavior data set
data2plot_behavior = data.behavior_bins{1,1};

% Events
% CS indexes
cs_trial_4hz = round(data.events{3, 1}./parameters.downsampling); 
cs_trial_8hz = round(data.events{4, 1}./parameters.downsampling); 

% Freezing indexes
freezing_start = data.behavior_bins{2, 1}(1,:);
freezing_end   = data.behavior_bins{2, 1}(1,:)+data.behavior_bins{2, 1}(2,:)-1;
% Time vector
behav_bins_time_v = linspace(0,length(data.behavior_bins{1,1})/(1/parameters.time_behav_bins),length(data.behavior_bins{1,1}));


py_path = "~/anaconda3/envs/Python_3_10/bin/python";
% colorbar off
% set(get(c,'ylabel'),'string','\fontsize{12} Zscore','Rotation',270);
% set(c,'XTickLabel',{'0',' ',' ',' ',' ','5'});
% view(0,90)

% Total power
normal_Spec = sum(abs(short_fft.data(closestfreq,:,:)));

title_ = {'mPFC IL';'mPFC PL';'dHPC'};

% Baseline power
%normal_Spec = abs(short_fft.data(closestfreq,(data.events{2, 1}(1,1)./parameters.decimated_srate)-21 : (data.events{2, 1}(1,1)./parameters.decimated_srate)-11,:));

f1 = figure;%('WindowState','maximized');
set(gcf,'color','w');
sc = [1,1,1800,880];
set(gcf, 'Position', sc);

for ii = 1:length(ch)
    
    subplot(2,3,ii)
    sgtitle({'Amplitude Spectrum via short-window FFT';['Hamming window = ' num2str(short_fft.timewin./1000) 's' ' - ' 'overlap = ' num2str(short_fft.overlap) '%']}) 
    contourf(short_fft.time,short_fft.freq(closestfreq),abs(short_fft.data(closestfreq,:,ch(ii)))./normal_Spec(:,:,ch(ii)),80,'linecolor','none');
    
    title(['channel ',num2str(ii+1)]);
    title(title_{ii})

    a = gca;
    a.TitleHorizontalAlignment = 'left';

    xlabel('Time (s)','FontSize',14), ylabel('Frequency (Hz)','FontSize',14)
    %xlim([short_fft.time(1) short_fft.time(end)])
    xlim([(data.events{3, 1}(1,1)./parameters.decimated_srate) - 20 (data.events{3, 1}(end,2)./parameters.decimated_srate) + 20])

    clim([-.007 .015])
    c = colorbar;
    c.Label.String = 'Normalized Power (A.U.)';
    Py_map = getPyPlot_cMap('RdGy_r', [], [], py_path);
    colormap(Py_map)

% % Events
%     z  = max(max(abs(short_fft.data(closestfreq_events,:,ch(ii)))));
%     zp = repmat (z, 1, length(short_fft.freq(closestfreq_events)));
%  
%     for jj = 1:length(short_fft.time_idx_t)
%         hold on
%         plot3(repmat(short_fft.time_idx_t(jj,1):short_fft.time_idx_t(jj,2),length(zp),1),freq2plot_events, zp,'w-','linew',2,'Color',[1, .4, .4,.6])
%     
%     end

        yyaxis right
        a = gca; % Get axis
        a.YColor = 'w';
        a.YTick = [];
        plot([data.timev_decimated(cs_trial_4hz(:,1));data.timev_decimated(cs_trial_4hz(:,2))], -50.*([ones(1,length(cs_trial_4hz(:,1))).*.5;ones(1,length(cs_trial_4hz(:,2))).*.5]),'k-','linew', 20,'Color',[.5, .7, 1])  
        hold all
        plot([behav_bins_time_v(freezing_start);behav_bins_time_v(freezing_end)], [ones(1,length(freezing_start)).*-450;ones(1,length(freezing_end)).*-450],'r-','linew', 2)
        ylim([-500 0])


end


for ii = 1:length(ch)

    subplot(2,3,ii+3)
    sgtitle({'Amplitude Spectrum via short-window FFT';['Hamming window = ' num2str(short_fft.timewin./1000) 's' ' - ' 'overlap = ' num2str(short_fft.overlap) '%']}) 
    contourf(short_fft.time,short_fft.freq(closestfreq),abs(short_fft.data(closestfreq,:,ch(ii)))./normal_Spec(:,:,ch(ii)),80,'linecolor','none');
    
    title(['channel ',num2str(ii+1)]);
    title(title_{ii})

    a = gca;
    a.TitleHorizontalAlignment = 'left';

    xlabel('Time (s)','FontSize',14), ylabel('Frequency (Hz)','FontSize',14)
    %xlim([short_fft.time(1) short_fft.time(end)])
    xlim([(data.events{4, 1}(1,1)./parameters.decimated_srate) - 20 (data.events{4, 1}(end,2)./parameters.decimated_srate) + 20])

    clim([-.007 .015])
    c = colorbar;
    c.Label.String = 'Normalized Power (A.U.)';
    Py_map = getPyPlot_cMap('RdGy_r', [], [], py_path);
    colormap(Py_map)

% % Events
%     z  = max(max(abs(short_fft.data(closestfreq_events,:,ch(ii)))));
%     zp = repmat (z, 1, length(short_fft.freq(closestfreq_events)));
%  
%     for jj = 1:length(short_fft.time_idx_t)
%         hold on
%         plot3(repmat(short_fft.time_idx_t(jj,1):short_fft.time_idx_t(jj,2),length(zp),1),freq2plot_events, zp,'w-','linew',2,'Color',[1, .4, .4,.6])
%     
%     end

        yyaxis right
        a = gca; % Get axis
        a.YColor = 'w';
        a.YTick = [];
        plot([data.timev_decimated(cs_trial_8hz(:,1));data.timev_decimated(cs_trial_8hz(:,2))], -50.*([ones(1,length(cs_trial_8hz(:,1))).*.5;ones(1,length(cs_trial_8hz(:,2))).*.5]),'k-','linew', 20,'Color',[0 0 .6])  
        hold all
        plot([behav_bins_time_v(freezing_start);behav_bins_time_v(freezing_end)], [ones(1,length(freezing_start)).*-450;ones(1,length(freezing_end)).*-450],'r-','linew', 3)
        ylim([-500 0])



end

clear ('f1','freq2plot_events','ch','steps','freq2plot','closestfreq','closestfreq_events','ii','z','zp','jj','a','c','normal_Spec')

%% Save

%newStr = regexprep(files.id.name,'.mat','_');
newStr1 = files.id(ms).name(1:end-12);
%path = files.FilesLoaded{1, 1}.folder;
path = '/Users/flavio/Desktop';

%name = strcat('E:\Projetos 2\Flavio\Samir\Analysis\Terceiro dia\',newStr1,newStr2,'_pw_mean_Trials_allCh');
name = strcat(path,'/',newStr1,'_short_fft');

%saveas(gcf,name,'png')

%set(gcf,'renderer','Painters')
exportgraphics(gcf,strcat(name,'.png'),'Resolution', 300)

%save(name,'short_fft','-v7.3')


close all

clear('name','newStr1','path')
%% last update 24/05/2024
%  listening: American Football : for sure

