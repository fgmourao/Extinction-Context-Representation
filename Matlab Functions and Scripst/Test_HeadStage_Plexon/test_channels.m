% Electrodos test


for dd = 1:size(data.test,1)

%Plot to check all channels in the same plot
data2plot_LFP = [];
%Choose data set
data2plot_LFP = data.test{dd, 1}.lfp{1, 1};

%Choose channels to plot
channels = 1:16;

%factor
factor = (channels)'*100;

%Set Figure
figure
set(gcf,'color','w');
sc = [1,1,1200,1200];
set(gcf, 'Position', sc);

box 'off'
hold on

subplot(8,4,[1 16])
%Select fields with data
r = plot(data.test{dd, 1}.timev, bsxfun(@plus, data2plot_LFP(channels,:), factor),'Color','[0.3, 0.3, 0.3]','linew',1);

a = gca; % Get axis
a.YColor = 'w';
a.YTick = [];
box off

a.XLim = [0 data.test{dd, 1}.timev(end)];

%xlabel('Time (s)')
ylabel('Channels 1 -> 16)')

title(files.id(dd).name)

%Clear trash
clear ('factor','channels','filter','substrate','session','str','sub','r','a','I','lh','lh_pos');

for ii = 17:32
    subplot(8,4,ii)
    plot(data.test{dd, 1}.timev,data2plot_LFP(ii-16,:),'Color','[0.3, 0.3, 0.3]','linew',1)
    xlim([0 data.test{dd, 1}.timev(end)]);
    ylim([-50 50]);
    title(['channel ',num2str(ii-16)]);

end

%% Save
newStr = files.id(dd).name(1:end-5);

path = '/Users/flavio/Desktop';
%path = files.FilesLoaded{1, 1}.folder;

name = strcat(path,'/',newStr,'_raw_plots');

% save data
saveas(gcf,name,'png') 
close all
clear('name','newStr','path') 

%% Spearman Correlation Matrix

[correlation.cormat_base{dd,1},correlation.PVAL_base{dd,1}] = corr(data2plot_LFP',data2plot_LFP','Type','s');

imagesc(correlation.cormat_base{dd,1});
xlabel ('channels')
ylabel ('channels')
xticks([1:16])
yticks([1:16])

title({files.id(dd).name;'Power Spearman`s correlation betwwen channels'})
colorbar
caxis([-1 1])

%% Save
newStr = files.id(dd).name(1:end-5);

path = '/Users/flavio/Desktop';
%path = files.FilesLoaded{1, 1}.folder;

name_1 = strcat(path,'/',newStr,'_correlation');

% save data
saveas(gcf,name_1,'png') 
close all
clear('name_1','newStr','path')


%%
% Initializing

% Time window
short_fft{dd,1}.timewin    = 25000; % in ms

% Convert time window to points
short_fft{dd,1}.timewinpnts  = hamming(round(short_fft{dd,1}.timewin/(1000/1000)));

% nFFT
%short_fft.nFFT = 2^nextpow2(round(short_fft.timewin/(1000/parameters.decimated_srate)));
short_fft{dd,1}.nFFT = 2^14; 

% Number of overlap samples
short_fft{dd,1}.overlap = 50;
short_fft{dd,1}.noverlap = floor(short_fft{dd,1}.overlap*0.01*(round(short_fft{dd,1}.timewin/(1000/1000))));


% Spectrogram
% lines: frequencies / columns: time / third dimension: channels

for ii = 1:size(data2plot_LFP,1)
    if ii == 1
       [short_fft{dd,1}.data(:,:,ii),short_fft{dd,1}.freq,short_fft{dd,1}.time] = spectrogram(data2plot_LFP(ii,:),short_fft{dd,1}.timewinpnts,short_fft{dd,1}.noverlap,short_fft{dd,1}.nFFT,1000);
    else
        short_fft{dd,1}.data(:,:,ii) = spectrogram(data2plot_LFP(ii,:),short_fft{dd,1}.timewinpnts,short_fft{dd,1}.noverlap,short_fft{dd,1}.nFFT,1000);
    end
end

clear ('ii','jj','not')


%% Plot to check full session. Channels per substrate 

% Choose channel
ch = 1:16;

%Define frequencies to plot in each subplot
steps = diff(short_fft{dd,1}.freq); % according to the fft time window

% For data
freq2plot = 1:steps(1):100;
closestfreq = dsearchn(short_fft{dd,1}.freq,freq2plot');


% Total power
%normal_Spec = sum(abs(short_fft.data{1, 1}(closestfreq,:,:)));


figure
set(gcf,'color','w');
sc = [1,1,1200,1200];
set(gcf, 'Position', sc);

for ii = 1:16
    subplot(4,4,ii)
%     Old function to add titles in subplot.  deprecated since 2021a
%     suptitle({'Amplitude Spectrum via short-window FFT';['(window = ' num2str(short_fft.timewin./1000) 's' ' - ' 'overlap = ' num2str(short_fft.overlap) '%)']}) 
% 
%     Add titles in subplot introduced in 18b
    sgtitle({'Amplitude Spectrum via short-window FFT';['Hamming window = ' num2str(short_fft{dd}.timewin./1000) 's' ' - ' 'overlap = ' num2str(short_fft{dd}.overlap) '%']}) 
    
    set(gcf,'color','white')

    contourf(short_fft{dd,1}.time,short_fft{dd,1}.freq(closestfreq),10*log10(abs(short_fft{dd,1}.data(closestfreq,:,ii))),80,'linecolor','none');
    
    title(['channel ',num2str(ii)]);
    a = gca;
    a.TitleHorizontalAlignment = 'left';

    xlabel('Time (s)','FontSize',14), ylabel('Frequency (Hz)','FontSize',14)
    xlim([short_fft{dd,1}.time(1) short_fft{dd,1}.time(end)])

    %clim([0 .10])
    c = colorbar;
    c.Label.String = 'Power (dB)';
    %c.Ticks =[0 .1];
   
%     if ii == 1
%         caxis([0 12*11^4])
%     elseif ii == 2
%         caxis([0 12*11^4])
%     elseif ii == 3
%         caxis([0 12*11^4])
%     elseif ii == 4
%         caxis([0 10*9^4])
%  

end

clear ('sc', 'f1','freq2plot_events','ch','steps','freq2plot','closestfreq','closestfreq_events','ii','z','zp','jj','a','c','normal_Spec')

%% Save
newStr = files.id(dd).name(1:end-5);

path = '/Users/flavio/Desktop';
%path = files.FilesLoaded{1, 1}.folder;

name_2 = strcat(path,'/',newStr,'_power_spectrum');

% save data
saveas(gcf,name_2,'png') 
close all
clear('name_2','newStr','path')

end
