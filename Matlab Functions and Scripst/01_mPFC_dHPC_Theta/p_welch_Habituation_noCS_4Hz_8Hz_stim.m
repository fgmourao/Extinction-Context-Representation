%%  Welch power spectral density estimate

% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  12/2023
% Last update: 01/2024


% First extract the data with: Extracting_LFPs_and_events.m
% ... then organize the data with the script: Pre_processing.m


%%
fprintf('\n Welch power spectral density estimate ... \n');

%% pwelch
%  [pxx,f] = pwelch(x,window,noverlap,f,fs)

% - cell
% - 1st column    - > baseline
% - 2nd column    - > cs events

% in each cell
% - rows          - > Hz
% - columns       - > channels
% 3th dimentions  - > CS-Trials

pw = [];

% Time window
pw.full_trial.parameters.baseline_timewin    = 2048; % in ms
pw.full_trial.parameters.STIM_Trials_timewin = 2048; % in ms
pw.full_trial.parameters.ITI_timewin         = 2048; % in ms

% Convert time window to points
pw.full_trial.parameters.baseline_timewinpnts     = hamming(round(pw.full_trial.parameters.baseline_timewin/(1000/parameters.decimated_srate)));
pw.full_trial.parameters.STIM_Trials_timewinpnts  = hamming(round(pw.full_trial.parameters.STIM_Trials_timewin/(1000/parameters.decimated_srate)));
pw.full_trial.parameters.ITI_timewinpnts          = hamming(round(pw.full_trial.parameters.ITI_timewin/(1000/parameters.decimated_srate)));

% nFFT
pw.full_trial.parameters.nFFT = 2^15; %4096; %2^nextpow2(pw.full_trial.parameters.baseline_timewinpnts));

% Number of overlap samples
pw.full_trial.parameters.overlap = 90;
pw.full_trial.parameters.baseline_noverlap = floor(pw.full_trial.parameters.overlap*0.01 * pw.full_trial.parameters.baseline_timewin);
pw.full_trial.parameters.STIM_Trials_noverlap = floor(pw.full_trial.parameters.overlap*0.01 * pw.full_trial.parameters.STIM_Trials_timewin);
pw.full_trial.parameters.ITI_noverlap = floor(pw.full_trial.parameters.overlap*0.01 * pw.full_trial.parameters.ITI_timewin);


% Baseline
not1 = 6;

for ii = 1:size(data.lfp{not1, 1},1)

    if ii == 1
        [pw.full_trial.Pxx{1,1}(ii,:),pw.full_trial.freq_baseline] = pwelch(data.lfp{not1, 1}(ii,:),pw.full_trial.parameters.baseline_timewinpnts,pw.full_trial.parameters.baseline_noverlap,pw.full_trial.parameters.nFFT,parameters.decimated_srate);

    else
        pw.full_trial.Pxx{1,1}(ii,:) = pwelch(data.lfp{not1, 1}(ii,:),pw.full_trial.parameters.baseline_timewinpnts,pw.full_trial.parameters.baseline_noverlap,pw.full_trial.parameters.nFFT,parameters.decimated_srate);

    end

end

clear ('ii')


% 4Hz or 8Hz opto stimulation
% Choose data from data.lfp according pre_processing.m define
not2 = 7;

for jj = 1:size(data.lfp{not2, 1},3)
    for ii = 1:size(data.lfp{not2, 1},1)
        if ii == 1
            [pw.full_trial.Pxx{2,1}(ii,:,jj),pw.full_trial.freq_STIM_Trials] = pwelch(data.lfp{not2, 1}(ii,:,jj),pw.full_trial.parameters.STIM_Trials_timewinpnts,pw.full_trial.parameters.STIM_Trials_noverlap,pw.full_trial.parameters.nFFT,parameters.decimated_srate);

        else
            pw.full_trial.Pxx{2,1}(ii,:,jj) = pwelch(data.lfp{not2, 1}(ii,:,jj),pw.full_trial.parameters.STIM_Trials_timewinpnts,pw.full_trial.parameters.STIM_Trials_noverlap,pw.full_trial.parameters.nFFT,parameters.decimated_srate);

        end

    end
end

clear ('jj','ii')

% % ITI or non-freezing epochs
% Choose data from data.lfp according pre_processing.m define
not3 = 8;

for jj = 1:size(data.lfp{not3, 1},3)
    for ii = 1:size(data.lfp{not3, 1},1)
        if ii == 1
            [pw.full_trial.Pxx{3,1}(ii,:,jj),pw.full_trial.freq_ITI] = pwelch(data.lfp{not3, 1}(ii,:,jj),pw.full_trial.parameters.ITI_timewinpnts,pw.full_trial.parameters.ITI_noverlap,pw.full_trial.parameters.nFFT,parameters.decimated_srate);

        else
            pw.full_trial.Pxx{3,1}(ii,:,jj) = pwelch(data.lfp{not3, 1}(ii,:,jj),pw.full_trial.parameters.ITI_timewinpnts,pw.full_trial.parameters.ITI_noverlap,pw.full_trial.parameters.nFFT,parameters.decimated_srate);

        end

    end
end

% 4Hz or 8Hz opto stimulation
% Choose data from data.lfp according pre_processing.m define
not4 = 9;

for jj = 1:size(data.lfp{not4, 1},3)
    for ii = 1:size(data.lfp{not4, 1},1)
        if ii == 1
            [pw.full_trial.Pxx{4,1}(ii,:,jj),pw.full_trial.freq_STIM_Trials] = pwelch(data.lfp{not4, 1}(ii,:,jj),pw.full_trial.parameters.STIM_Trials_timewinpnts,pw.full_trial.parameters.STIM_Trials_noverlap,pw.full_trial.parameters.nFFT,parameters.decimated_srate);

        else
            pw.full_trial.Pxx{4,1}(ii,:,jj) = pwelch(data.lfp{not4, 1}(ii,:,jj),pw.full_trial.parameters.STIM_Trials_timewinpnts,pw.full_trial.parameters.STIM_Trials_noverlap,pw.full_trial.parameters.nFFT,parameters.decimated_srate);

        end

    end
end

clear ('jj','ii')

% % ITI or non-freezing epochs
% Choose data from data.lfp according pre_processing.m define
not5 = 10;

for jj = 1:size(data.lfp{not5, 1},3)
    for ii = 1:size(data.lfp{not5, 1},1)
        if ii == 1
            [pw.full_trial.Pxx{5,1}(ii,:,jj),pw.full_trial.freq_ITI] = pwelch(data.lfp{not5, 1}(ii,:,jj),pw.full_trial.parameters.ITI_timewinpnts,pw.full_trial.parameters.ITI_noverlap,pw.full_trial.parameters.nFFT,parameters.decimated_srate);

        else
            pw.full_trial.Pxx{5,1}(ii,:,jj) = pwelch(data.lfp{not5, 1}(ii,:,jj),pw.full_trial.parameters.ITI_timewinpnts,pw.full_trial.parameters.ITI_noverlap,pw.full_trial.parameters.nFFT,parameters.decimated_srate);

        end

    end
end

clear ('jj','ii')

% figure
% plot(pw.full_trial.freq_baseline,pw.full_trial.Pxx{1,1}(1,:))
% hold on
% plot(pw.full_trial.freq_STIM_Trials,mean(pw.full_trial.Pxx{2,1}(1,:,1:3),3))
% plot(pw.full_trial.freq_STIM_Trials,mean(pw.full_trial.Pxx{3,1}(1,:,1:3),3))
%
% xlim([3 12])

clear ('not1','not2','not3')

%% Save

%newStr = regexprep(files.id.name,'.mat','_');
%newStr = files.id(ms).name(1:end-8);
newStr = id(1:end-5);

path = '/Users/flavio/Desktop';
%path = files.FilesLoaded{1, 1}.folder;

name = strcat(path,'/',newStr,'_pw_Total_Power');

% save data
save(name,'pw','-v7.3')

clear('name','newStr','path')

%% Set theta range

steps                                    = diff(pw.full_trial.freq_baseline); % according to the fft time window

pw.full_trial.parameters.frex_2_12Hz     = 2:steps(1):12;
pw.full_trial.parameters.frex_idx_2_12Hz = dsearchn(pw.full_trial.freq_baseline,pw.full_trial.parameters.frex_2_12Hz');

clear('steps')

%% Plot to check - normalize by Total Power

% choose channel
ch1 = 1;
ch2 = 2;
ch3 = 3;

f = figure;
set(gcf, 'Position', get(0, 'Screensize'));
set(gcf,'color','w');


% Baseline
b(1) = subplot(3,7,1);
plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),...
    pw.full_trial.Pxx{1, 1}(ch1,pw.full_trial.parameters.frex_idx_2_12Hz)./sum(pw.full_trial.Pxx{1, 1}(ch1,pw.full_trial.parameters.frex_idx_2_12Hz),2),'linewidth',2,'color',[.6 .6 .6] )
xlim([2 12])
ylabel('Amplitude')
xlabel('Hz')
box off
title('baseline')


b(2) = subplot(3,7,8);
plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),...
    pw.full_trial.Pxx{1, 1}(ch2,pw.full_trial.parameters.frex_idx_2_12Hz)./sum(pw.full_trial.Pxx{1, 1}(ch2,pw.full_trial.parameters.frex_idx_2_12Hz),2),'linewidth',2,'color',[.6 .6 .6] )
xlim([2 12])
ylabel('Amplitude')
box off
xlabel('Hz')

b(3) = subplot(3,7,15);
plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),...
    pw.full_trial.Pxx{1, 1}(ch3,pw.full_trial.parameters.frex_idx_2_12Hz)./sum(pw.full_trial.Pxx{1, 1}(ch3,pw.full_trial.parameters.frex_idx_2_12Hz),2),'linewidth',2,'color',[.6 .6 .6] )
xlim([2 12])
box off
ylabel('Amplitude')
xlabel('Hz')



% STM Block 1

for ii = 1:size(pw.full_trial.Pxx{2, 1},3)
    t(ii) = subplot(3,7,ii+1);

    plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),...
        pw.full_trial.Pxx{2, 1}(ch1,pw.full_trial.parameters.frex_idx_2_12Hz,ii)./sum(pw.full_trial.Pxx{2, 1}(ch1,pw.full_trial.parameters.frex_idx_2_12Hz),2),'linewidth',2,'Color',[.5, .7, 1])
    hold all
    plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),...
        pw.full_trial.Pxx{4, 1}(ch1,pw.full_trial.parameters.frex_idx_2_12Hz,ii)./sum(pw.full_trial.Pxx{4, 1}(ch1,pw.full_trial.parameters.frex_idx_2_12Hz),2),'linewidth',2,'Color',[0 0 .6])

    xlim([2 12])
    ylabel('Amplitude')
    xlabel('Hz')
    box off
    title(['STIM ' num2str(ii)])


    t(ii+8) = subplot(3,7,ii+8);
    plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),...
        pw.full_trial.Pxx{2, 1}(ch2,pw.full_trial.parameters.frex_idx_2_12Hz,ii)./sum(pw.full_trial.Pxx{2, 1}(ch2,pw.full_trial.parameters.frex_idx_2_12Hz),2),'linewidth',2,'Color',[.5, .7, 1])
    hold all
    plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),...
        pw.full_trial.Pxx{4, 1}(ch2,pw.full_trial.parameters.frex_idx_2_12Hz,ii)./sum(pw.full_trial.Pxx{4, 1}(ch2,pw.full_trial.parameters.frex_idx_2_12Hz),2),'linewidth',2,'Color',[0 0 .6])

    xlim([2 12])
    ylabel('Amplitude')
    box off
    xlabel('Hz')

    t(ii+15) = subplot(3,7,ii+15);
    plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),...
        pw.full_trial.Pxx{2, 1}(ch3,pw.full_trial.parameters.frex_idx_2_12Hz,ii)./sum(pw.full_trial.Pxx{2, 1}(ch3,pw.full_trial.parameters.frex_idx_2_12Hz),2),'linewidth',2,'Color',[.5, .7, 1])
    hold all
    plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),...
        pw.full_trial.Pxx{4, 1}(ch3,pw.full_trial.parameters.frex_idx_2_12Hz,ii)./sum(pw.full_trial.Pxx{4, 1}(ch3,pw.full_trial.parameters.frex_idx_2_12Hz),2),'linewidth',2,'Color',[0 0 .6])

    xlim([2 12])
    box off
    ylabel('Amplitude')
    xlabel('Hz')

end

% STM Block mean - ONLY BEST TRIALS 

m(1) = subplot(3,7,7);
plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),...
    mean(pw.full_trial.Pxx{2, 1}(ch1,pw.full_trial.parameters.frex_idx_2_12Hz,CSIT{ms}(1,:))./sum(pw.full_trial.Pxx{2, 1}(ch1,pw.full_trial.parameters.frex_idx_2_12Hz,CSIT{ms}(1,:)),2),3),'linewidth',2,'Color',[.5, .7, 1])
hold all
plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),...
    mean(pw.full_trial.Pxx{4, 1}(ch1,pw.full_trial.parameters.frex_idx_2_12Hz,CSIT{ms}(2,:))./sum(pw.full_trial.Pxx{4, 1}(ch1,pw.full_trial.parameters.frex_idx_2_12Hz,CSIT{ms}(2,:)),2),3),'linewidth',2,'Color',[0 0 .6])

xlim([2 12])
ylabel('Amplitude')
xlabel('Hz')
box off
title(['Averaged Trials'])


m(2) = subplot(3,7,14);
plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),...
    mean(pw.full_trial.Pxx{2, 1}(ch2,pw.full_trial.parameters.frex_idx_2_12Hz,CSIT{ms}(1,:))./sum(pw.full_trial.Pxx{2, 1}(ch2,pw.full_trial.parameters.frex_idx_2_12Hz,CSIT{ms}(1,:)),2),3),'linewidth',2,'Color',[.5, .7, 1])
hold all
plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),...
    mean(pw.full_trial.Pxx{4, 1}(ch2,pw.full_trial.parameters.frex_idx_2_12Hz,CSIT{ms}(2,:))./sum(pw.full_trial.Pxx{4, 1}(ch2,pw.full_trial.parameters.frex_idx_2_12Hz,CSIT{ms}(2,:)),2),3),'linewidth',2,'Color',[0 0 .6])

xlim([2 12])
ylabel('Amplitude')
box off
xlabel('Hz')

m(3) = subplot(3,7,21);
plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),...
    mean(pw.full_trial.Pxx{2, 1}(ch3,pw.full_trial.parameters.frex_idx_2_12Hz,CSIT{ms}(1,:))./sum(pw.full_trial.Pxx{2, 1}(ch3,pw.full_trial.parameters.frex_idx_2_12Hz,CSIT{ms}(1,:)),2),3),'linewidth',2,'Color',[.5, .7, 1])
hold all
plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),...
    mean(pw.full_trial.Pxx{4, 1}(ch3,pw.full_trial.parameters.frex_idx_2_12Hz,CSIT{ms}(2,:))./sum(pw.full_trial.Pxx{4, 1}(ch3,pw.full_trial.parameters.frex_idx_2_12Hz,CSIT{ms}(2,:)),2),3),'linewidth',2,'Color',[0 0 .6])

xlim([2 12])
box off
ylabel('Amplitude')
xlabel('Hz')


linkaxes([b t m],'xy');
b(1).YLim  = [0 0.045];

clear ('ii','ch1','ch2','ch3')

%% Save

newStr = id(1:end-5);

path = '/Users/flavio/Desktop';
%path = files.FilesLoaded{1, 1}.folder;

name = strcat(path,'/',newStr,'_pw_Total_Power_best_Channels_averaged');

saveas(gcf,name,'png')

set(gcf,'renderer','Painters')
exportgraphics(gcf,strcat(name,'.eps'),'Resolution', 300)

close all

clear('name','newStr1','path')

%% Plot to check - baseline norm

% choose channel
ch1 = 1;
ch2 = 2;
ch3 = 3;

figure
set(gcf,'color','w');
sc = [1,1,1800,580];
set(gcf, 'Position', sc);


% STM Block Zscore - ONLY BEST TRIALS 

bn(1) = subplot(1,3,1);
plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),...
        mean((pw.full_trial.Pxx{2, 1}(ch1,pw.full_trial.parameters.frex_idx_2_12Hz,CSIT{ms}(1,:)) - pw.full_trial.Pxx{1, 1}(ch1,pw.full_trial.parameters.frex_idx_2_12Hz))...
    ./pw.full_trial.Pxx{1, 1}(ch1,pw.full_trial.parameters.frex_idx_2_12Hz),3),'linewidth',2,'color',[.5, .7, 1] )
hold all
plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),...
        mean((pw.full_trial.Pxx{4, 1}(ch1,pw.full_trial.parameters.frex_idx_2_12Hz,CSIT{ms}(2,:)) - pw.full_trial.Pxx{1, 1}(ch1,pw.full_trial.parameters.frex_idx_2_12Hz))...
    ./pw.full_trial.Pxx{1, 1}(ch1,pw.full_trial.parameters.frex_idx_2_12Hz),3),'linewidth',2,'color',[0, 0, .6] )

title('mPFC IL','FontSize',16)
xlim([2 12])
%ylim([-2 5])

ax=gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;

ylabel({'Amplitude Spectrum';'(Change from baseline)'},'FontSize',16)
xlabel('Hz','FontSize',16)
box off


bn(2) = subplot(1,3,2);
plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),...
        mean((pw.full_trial.Pxx{2, 1}(ch2,pw.full_trial.parameters.frex_idx_2_12Hz,CSIT{ms}(1,:)) - pw.full_trial.Pxx{1, 1}(ch2,pw.full_trial.parameters.frex_idx_2_12Hz))...
    ./pw.full_trial.Pxx{1, 1}(ch2,pw.full_trial.parameters.frex_idx_2_12Hz),3),'linewidth',2,'color',[.5, .7, 1] )
hold all
plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),...
        mean((pw.full_trial.Pxx{4, 1}(ch2,pw.full_trial.parameters.frex_idx_2_12Hz,CSIT{ms}(2,:)) - pw.full_trial.Pxx{1, 1}(ch2,pw.full_trial.parameters.frex_idx_2_12Hz))...
    ./pw.full_trial.Pxx{1, 1}(ch2,pw.full_trial.parameters.frex_idx_2_12Hz),3),'linewidth',2,'color',[0, 0, .6] )

title('mPFC PL','FontSize',16)
xlim([2 12])
%ylim([-2 5])

ax=gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;

box off
xlabel('Hz','FontSize',16)


bn(3) = subplot(1,3,3);
plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),...
        mean((pw.full_trial.Pxx{2, 1}(ch3,pw.full_trial.parameters.frex_idx_2_12Hz,CSIT{ms}(1,:)) - pw.full_trial.Pxx{1, 1}(ch3,pw.full_trial.parameters.frex_idx_2_12Hz))...
    ./pw.full_trial.Pxx{1, 1}(ch3,pw.full_trial.parameters.frex_idx_2_12Hz),3),'linewidth',2,'color',[.5, .7, 1] )
hold all
plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),...
        mean((pw.full_trial.Pxx{4, 1}(ch3,pw.full_trial.parameters.frex_idx_2_12Hz,CSIT{ms}(2,:)) - pw.full_trial.Pxx{1, 1}(ch3,pw.full_trial.parameters.frex_idx_2_12Hz))...
    ./pw.full_trial.Pxx{1, 1}(ch3,pw.full_trial.parameters.frex_idx_2_12Hz),3),'linewidth',2,'color',[0, 0, .6] )

title('dHPC','FontSize',16)
xlim([2 12])
%ylim([-2 5])

ax=gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;

box off
xlabel('Hz','FontSize',16)

legend('4Hz','8Hz','FontSize',16)

linkaxes(bn,'xy');
bn(1).YLim  = [-1 8];

clear ('ii','ch1','ch2','ch3')

%% Save

newStr = id(1:end-5);

path = '/Users/flavio/Desktop';
%path = files.FilesLoaded{1, 1}.folder;

name = strcat(path,'/',newStr,'_pw_Total_Power_Baseline_norm');

saveas(gcf,name,'png')

set(gcf,'renderer','Painters')
exportgraphics(gcf,strcat(name,'.eps'),'Resolution', 300)

close all

clear('name','newStr1','path')
%% last update 20/05/2024
%  listening: Sonic Youth - Disconnection Notice
