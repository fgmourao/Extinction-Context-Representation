
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
short_fft.timewin    = 6000; % in ms

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

%% Plot to check full session. Channels per substrate 

% Choose channel
ch = 1:3;

%Define frequencies to plot in each subplot
steps = diff(short_fft.freq); % according to the fft time window

% For data
freq2plot = 2:steps(1):12;
closestfreq = dsearchn(short_fft.freq,freq2plot');

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

for ii = 1:3
    
    subplot(2,3,ii)
    

    sgtitle({'Amplitude Spectrum via short-window FFT';['Hamming window = ' num2str(short_fft.timewin./1000) 's' ' - ' 'overlap = ' num2str(short_fft.overlap) '%']}) 
    contourf(short_fft.time,short_fft.freq(closestfreq),abs(short_fft.data(closestfreq,:,ch(ii)))./normal_Spec(:,:,ch(ii)),80,'linecolor','none');
    
    title(['channel ',num2str(ii+1)]);
    title(title_{ii})

    a = gca;
    a.TitleHorizontalAlignment = 'left';

    xlabel('Time (s)','FontSize',14), ylabel('Frequency (Hz)','FontSize',14)
    %xlim([short_fft.time(1) short_fft.time(end)])

    clim([0 0.006])
    c = colorbar;
    c.Label.String = 'Normalized Power (A.U.)';
    Py_map = getPyPlot_cMap('RdBu_r', [], [], py_path);
    colormap(Py_map)

    xlim([193-3 203+3])


        yyaxis right
%         plot(behav_bins_time_v, movmean(data2plot_behavior,(1/parameters.time_behav_bins))-500,'linew', 1, 'Color',[0.6, 0.6, 0.6, 0.6])

        a = gca; % Get axis
%         a.YColor = 'w';
%         a.YTick = [];
        hold all
        plot([behav_bins_time_v(freezing_start);behav_bins_time_v(freezing_end)], [ones(1,length(freezing_start)).*-425;ones(1,length(freezing_end)).*-425],'k-','linew', 2,'Color',[.6, 0, 0])
        %ylim([-2000 0])
        xlim([180-3 190+3])

        xlim([180-3 190+3])

        subplot(2,3,ii+3)
        plot(behav_bins_time_v, data2plot_behavior,'linew', 1, 'Color',[0.6, 0.6, 0.6, 0.6])
        xlim([180-3 190+3])
        %xlim([data.events{2, 1}(1,1)/1000 data.events{2, 1}(1,2)/1000])
end




clear ('f1','freq2plot_events','ch','steps','freq2plot','closestfreq','closestfreq_events','z','zp','jj','a','c','normal_Spec')

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

