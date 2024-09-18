%%  Welch power spectral density estimate

% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  12/2023
% Last update: 05/2024


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
        [pw.full_trial.Pxx{1,1}(ii,:),pw.full_trial.freq_baseline] = pwelch(data.lfp{not1, 1}(ii,B_clean{ms}(1):B_clean{ms}(2)),pw.full_trial.parameters.baseline_timewinpnts,pw.full_trial.parameters.baseline_noverlap,pw.full_trial.parameters.nFFT,parameters.decimated_srate);

    else
        pw.full_trial.Pxx{1,1}(ii,:) = pwelch(data.lfp{not1, 1}(ii,B_clean{ms}(1):B_clean{ms}(2)),pw.full_trial.parameters.baseline_timewinpnts,pw.full_trial.parameters.baseline_noverlap,pw.full_trial.parameters.nFFT,parameters.decimated_srate);

    end

end

clear ('ii')


% CS-Tone
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


% % ITI
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


% figure
% plot(pw.full_trial.freq_baseline,pw.full_trial.Pxx{1,1}(1,:))
% hold on
% plot(pw.full_trial.freq_STIM_Trials,mean(pw.full_trial.Pxx{2,1}(1,:,1:3),3))
% plot(pw.full_trial.freq_STIM_Trials,mean(pw.full_trial.Pxx{3,1}(1,:,1:3),3))
%
% xlim([3 12])


%% Set theta range

steps                                    = diff(pw.full_trial.freq_baseline); % according to the fft time window

pw.full_trial.parameters.frex_2_12Hz     = 2:steps(1):12;
pw.full_trial.parameters.frex_idx_2_12Hz = dsearchn(pw.full_trial.freq_baseline,pw.full_trial.parameters.frex_2_12Hz');

clear('steps')


%% Normalization and trials average

pw.full_trial.Pxx_Baseline_norm         = [];
pw.full_trial.Pxx_TotalPower_norm       = [];
pw.full_trial.Pxx_Baseline_norm_mean    = [];
pw.full_trial.Pxx_TotalPower_norm_mean  = [];


for jj = 1:size(pw.full_trial.Pxx,1)

    pw.full_trial.Pxx_Baseline_norm{jj,1}   = pw.full_trial.Pxx{jj,1}(:,pw.full_trial.parameters.frex_idx_2_12Hz,:)./pw.full_trial.Pxx{1, 1}(:,pw.full_trial.parameters.frex_idx_2_12Hz); % Baseline  normalization
    pw.full_trial.Pxx_TotalPower_norm{jj,1} = pw.full_trial.Pxx{jj,1}(:,pw.full_trial.parameters.frex_idx_2_12Hz,:)./sum(pw.full_trial.Pxx{jj,1}(:,pw.full_trial.parameters.frex_idx_2_12Hz,:),2); % total power normalization


    if jj == 1 % baseline condition
        pw.full_trial.Pxx_Baseline_norm_mean{jj,1}   = pw.full_trial.Pxx_Baseline_norm{jj,1};   % baseline cell will in 1 of course....but just to keep the organization
        pw.full_trial.Pxx_TotalPower_norm_mean{jj,1} = pw.full_trial.Pxx_TotalPower_norm{jj,1}; % baseline cell will in 1 of course....but just to keep the organization

    else

        pw.full_trial.Pxx_Baseline_norm_mean{jj,1}    = mean(pw.full_trial.Pxx_Baseline_norm{jj,1},3);
        pw.full_trial.Pxx_TotalPower_norm_mean{jj,1}  = mean(pw.full_trial.Pxx_TotalPower_norm{jj,1},3);
            
    end

end


%% Plot to check - for now I am considring just CS-Tone + 8Hz opto stimulation (10s opto + CS-Tone)

%data_2_plot = pw.full_trial.Pxx_stats_baseline;
data_2_plot      = pw.full_trial.Pxx_TotalPower_norm;
data_2_plot_mean = pw.full_trial.Pxx_TotalPower_norm_mean;

% choose channel
ch1 = 1;
ch2 = 2;
ch3 = 3;

f = figure;
set(gcf, 'Position', get(0, 'Screensize'));
set(gcf,'color','w');


% Baseline
b(1) = subplot(3,12,1);
plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),data_2_plot{1, 1}(ch1,:),'linewidth',2,'color',[.6 .6 .6] )
xlim([2 12])
ylabel({'mPFC (IL)';'Amplitude'})
xlabel('Hz')
box off
title({'Baseline';[]})


b(2) = subplot(3,12,13);
plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),data_2_plot{1, 1}(ch2,:),'linewidth',2,'color',[.6 .6 .6] )
xlim([2 12])
ylabel({'mPFC (PL)';'Amplitude'})
box off
xlabel('Hz')

b(3) = subplot(3,12,25);
plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),data_2_plot{1, 1}(ch3,:),'linewidth',2,'color',[.6 .6 .6] )

xlim([2 12])
box off
ylabel({'dHPC';'Amplitude'})
xlabel('Hz')



% CS-Tones

for ii = 1: size(data_2_plot{2, 1},3)

    t(ii) = subplot(3,12,ii+1);
    plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),data_2_plot{2, 1}(ch1,:,ii),'linewidth',2,'Color',[0 0 .6])

    xlim([2 12])
    %ylabel('Amplitude')
    xlabel('Hz')
    box off
    title({['Stim ' num2str(ii)];[]})

    t(ii+11) = subplot(3,12,ii+13);
    plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),data_2_plot{2, 1}(ch2,:,ii),'linewidth',2,'Color',[0 0 .6])

    xlim([2 12])
    %ylabel('Amplitude')
    box off
    xlabel('Hz')

    t(ii+21) = subplot(3,12,ii+25);
    plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),data_2_plot{2, 1}(ch3,:,ii),'linewidth',2,'Color',[0 0 .6])

    xlim([2 12])
    box off
    %ylabel('Amplitude')
    xlabel('Hz')

end


% Trial averaged

m(1) = subplot(3,12,12);
plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),data_2_plot_mean{2, 1}(ch1,:),'linewidth',2,'color',[.1 .1 .1] )
xlim([2 12])
xlabel('Hz')
box off
title({'Averaged Trials';[]})


m(2) = subplot(3,12,24);
plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),data_2_plot_mean{2, 1}(ch2,:),'linewidth',2,'color',[.1 .1 .1] )
xlim([2 12])
box off
xlabel('Hz')

m(3) = subplot(3,12,36);
plot(pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12Hz),data_2_plot_mean{2, 1}(ch3,:),'linewidth',2,'color',[.1 .1 .1] )

xlim([2 12])
box off
xlabel('Hz')

linkaxes([b t m],'xy');
%b(1).YLim  = [0 15];
b(1).YLim  = [0 0.020];

clear ('ii','ch1','ch2','ch3')

%% Save figure

newStr = id(1:end-5);

path = '/Users/flavio/Desktop';
%path = files.FilesLoaded{1, 1}.folder;

%name = strcat(path,'/',newStr,'_pw_Total_Power_best_5_min_blocks_baseline_norm');
name = strcat(path,'/',newStr,'_pw_Total_Power_norm');

saveas(gcf,name,'png')

set(gcf,'renderer','Painters')
exportgraphics(gcf,strcat(name,'.eps'),'Resolution', 300)

close all

clear('name','newStr1','path')

%% Save data

%newStr = regexprep(files.id.name,'.mat','_');
%newStr = files.id(ms).name(1:end-8);
newStr = id(1:end-5);

path = '/Users/flavio/Desktop';
%path = files.FilesLoaded{1, 1}.folder;

name = strcat(path,'/',newStr,'_pw_Baseline_and_Total_Power');

% save data
save(name,'pw','-v7.3')

clear('name','newStr','path')
%% last update 29/05/2024
%  listening:
