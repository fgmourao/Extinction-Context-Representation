%%  Welch power spectral density estimate 

% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  12/2023
% Last update: 01/2024

%%
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
pw.full_trial.baseline_timewin    = 1000; % in ms
pw.full_trial.CS_Trials_timewin   = 1000; % in ms
pw.full_trial.ITI_timewin         = 1000; % in ms

% Convert time window to points
pw.full_trial.baseline_timewinpnts   = hamming(round(pw.full_trial.baseline_timewin/(1000/parameters.decimated_srate)));
pw.full_trial.CS_Trials_timewinpnts  = hamming(round(pw.full_trial.CS_Trials_timewin/(1000/parameters.decimated_srate)));
pw.full_trial.ITI_timewinpnts        = hamming(round(pw.full_trial.ITI_timewin/(1000/parameters.decimated_srate)));

% nFFT
pw.full_trial.nFFT = 2^15;%4096; %2^nextpow2(pw.full_trial.baseline_timewinpnts));

% Number of overlap samples
pw.full_trial.overlap = 50;
pw.full_trial.baseline_noverlap = floor(pw.full_trial.overlap*0.01 * pw.full_trial.baseline_timewin);
pw.full_trial.CS_Trials_noverlap = floor(pw.full_trial.overlap*0.01 * pw.full_trial.CS_Trials_timewin);
pw.full_trial.ITI_noverlap = floor(pw.full_trial.overlap*0.01 * pw.full_trial.ITI_timewin);

% Baseline
% Choose data from data.lfp according pre_processing.m
not1 = 6;

for ii = 1:size(data.lfp{not1, 1},1)

    if ii == 1
        [pw.full_trial.Pxx{1,1}(ii,:),pw.full_trial.freq_baseline] = pwelch(data.lfp{not1, 2}(ii,:),pw.full_trial.baseline_timewinpnts,pw.full_trial.baseline_noverlap,pw.full_trial.nFFT,parameters.decimated_srate);

    else
        pw.full_trial.Pxx{1,1}(ii,:) = pwelch(data.lfp{not1, 1}(ii,:),pw.full_trial.baseline_timewinpnts,pw.full_trial.baseline_noverlap,pw.full_trial.nFFT,parameters.decimated_srate);

    end

end

clear ('ii')


% CS-Trials
% Choose data from data.lfp according pre_processing.m define
not2 = 7;

for jj = 1:size(data.lfp{not2, 1},3)
    for ii = 1:size(data.lfp{not2, 1},1)
        if ii == 1
            [pw.full_trial.Pxx{2,1}(ii,:,jj),pw.full_trial.freq_CS_trials] = pwelch(data.lfp{not2, 2}(ii,:,jj),pw.full_trial.CS_Trials_timewinpnts,pw.full_trial.CS_Trials_noverlap,pw.full_trial.nFFT,parameters.decimated_srate);

        else
            pw.full_trial.Pxx{2,1}(ii,:,jj) = pwelch(data.lfp{not2, 1}(ii,:,jj),pw.full_trial.CS_Trials_timewinpnts,pw.full_trial.CS_Trials_noverlap,pw.full_trial.nFFT,parameters.decimated_srate);

        end

    end
end

clear ('jj','ii')

% ITI
% Choose data from data.lfp according pre_processing.m define
not3 = 8;

for jj = 1:size(data.lfp{not3, 1},3)
    for ii = 1:size(data.lfp{not3, 1},1)
        if ii == 1
            [pw.full_trial.Pxx{3,1}(ii,:,jj),pw.full_trial.freq_ITI] = pwelch(data.lfp{not3, 2}(ii,:,jj),pw.full_trial.ITI_timewinpnts,pw.full_trial.ITI_noverlap,pw.full_trial.nFFT,parameters.decimated_srate);

        else
            pw.full_trial.Pxx{3,1}(ii,:,jj) = pwelch(data.lfp{not3, 1}(ii,:,jj),pw.full_trial.ITI_timewinpnts,pw.full_trial.ITI_noverlap,pw.full_trial.nFFT,parameters.decimated_srate);

        end

    end
end

clear ('jj','ii')

% figure
% plot(pw.full_trial.freq_baseline,pw.full_trial.Pxx{1,1}(1,:))
% hold on
% plot(pw.full_trial.freq_CS_trials,mean(pw.full_trial.Pxx{2,1}(1,:,1:3),3))
% plot(pw.full_trial.freq_CS_trials,mean(pw.full_trial.Pxx{3,1}(1,:,1:3),3))
% 
% xlim([3 12])

clear ('not1','not2','not3')

%% Extract data
% Baseline Normalization 

% Define band frequency
% 2 - 12 Hertz
steps         = diff(pw.full_trial.freq_baseline); % according to the fft time window
freq2plot     = 2:steps(1):12;
pw.stats_baseline{5,1}  = dsearchn(pw.full_trial.freq_baseline,freq2plot');

% CS-Trials
pw.stats_baseline{1,1} = pw.full_trial.Pxx{2,1}(:,pw.stats_baseline{5,1},:)./pw.full_trial.Pxx{1,1}(:,pw.stats_baseline{5,1}); % / Baseline
pw.stats_baseline{2,1} = 10*log10(pw.full_trial.Pxx{2,1}(:,pw.stats_baseline{5,1},:)./pw.full_trial.Pxx{1,1}(:,pw.stats_baseline{5,1})); % log power / Baseline
pw.stats_baseline{3,1} = (pw.full_trial.Pxx{2,1}(:,pw.stats_baseline{5,1},:) - pw.full_trial.Pxx{1,1}(:,pw.stats_baseline{5,1}))./pw.full_trial.Pxx{1,1}(:,pw.stats_baseline{5,1}); % / percentage
pw.stats_baseline{4,1} = (pw.full_trial.Pxx{2,1}(:,pw.stats_baseline{5,1},:) - mean(pw.full_trial.Pxx{1,1}(:,pw.stats_baseline{5,1}),2))./ std(pw.full_trial.Pxx{1,1}(:,pw.stats_baseline{5,1}),[],2); % / Z-score

clear('steps','freq2plot')

%% Extract data
% Total Power Normalization according to the Maren & Totty paper

% First line: Full spectrum between 1 and 12 Hz
% Second line: "3 to 6 hertz"
% Third line: "6 to 9 hertz"

% Fifth line: "3 to 6 hertz" -> Data related to Extinction. Following the laboratory method, averages were calculated for every 5 trials.
% sixth line: "6 to 9 hertz" -> Data related to Extinction. Following the laboratory method, averages were calculated for every 5 trials.


% First column: "Baseline"
% Second column: "CS-Trials"
% Third column: "ITI"

% Within each cell:
%  - Rows: Channels
%  - Columns: Trials

% Define band frequency
% 2 - 12 Hertz / Keep fullspec to plot
steps         = diff(pw.full_trial.freq_baseline); % according to the fft time window
freq2plot_1   = 2:steps(1):12;
closestfreq_1 = dsearchn(pw.full_trial.freq_baseline,freq2plot_1');

pw.stats_total_power{1,1} = pw.full_trial.Pxx{1,1}(:,closestfreq_1)./sum(pw.full_trial.Pxx{1,1}(:,closestfreq_1),2); % Baseline
pw.stats_total_power{1,2} = pw.full_trial.Pxx{2,1}(:,closestfreq_1,:)./sum(pw.full_trial.Pxx{2,1}(:,closestfreq_1),2); % CS-Trials
pw.stats_total_power{1,3} = pw.full_trial.Pxx{3,1}(:,closestfreq_1,:)./sum(pw.full_trial.Pxx{3,1}(:,closestfreq_1),2); % ITI


% Define band frequency
% 3 - 6 Hertz / average
steps         = diff(freq2plot_1); % according to the fft time window
freq2plot_2   = 3:steps(1):6;
closestfreq_2 = dsearchn(freq2plot_1',freq2plot_2');

pw.stats_total_power{2,1} = mean(pw.stats_total_power{1,1}(:,closestfreq_2),2); % Baseline
pw.stats_total_power{2,2} = squeeze(mean(pw.stats_total_power{1,2}(:,closestfreq_2,:),2)); % CS-Trials
pw.stats_total_power{2,3} = squeeze(mean(pw.stats_total_power{1,3}(:,closestfreq_2,:),2)); % ITI

% Graph based on Maren`s paper - nature communication
% Data related to Extinction. Following the laboratory method, averages were calculated for every 5 trials.
% pw.stats_total_power{5,2} = squeeze(mean(reshape(pw.stats_total_power{2,2},size(pw.stats_total_power{2,2},1),[],9),2)); % CS-Trials
% pw.stats_total_power{5,3} = squeeze(mean(reshape(pw.stats_total_power{2,3},size(pw.stats_total_power{2,3},1),[],9),2)); % ITI


% Define band frequency
% 6 - 9 Hertz / average
steps         = diff(pw.full_trial.freq_baseline); % according to the fft time window
freq2plot_3   = 6:steps(1):9;
closestfreq_3 = dsearchn(freq2plot_1',freq2plot_3');

pw.stats_total_power{3,1} = mean(pw.stats_total_power{1,1}(:,closestfreq_3),2); % Baseline
pw.stats_total_power{3,2} = squeeze(mean(pw.stats_total_power{1,2}(:,closestfreq_3,:),2)); % CS-Trials
pw.stats_total_power{3,3} = squeeze(mean(pw.stats_total_power{1,3}(:,closestfreq_3,:),2)); % ITI

% Graph based on Maren`s paper - nature communication
% Data related to Extinction. Following the laboratory method, averages were calculated for every 5 trials.
% pw.stats_total_power{6,2} = squeeze(mean(reshape(pw.stats_total_power{3,2},size(pw.stats_total_power{3,2},1),[],9),2)); % CS-Trials
% pw.stats_total_power{6,3} = squeeze(mean(reshape(pw.stats_total_power{3,3},size(pw.stats_total_power{3,3},1),[],9),2)); % ITI


clear('steps','freq2plot_1','closestfreq_1','freq2plot_2','closestfreq_2','freq2plot_3','closestfreq_3')

%% Plot to check 5 block session average. Averaging Channels per substrate as well
% Total Power normalization

data_2_plot_1 = pw.stats_total_power{1,1};
data_2_plot_2 = pw.stats_total_power{1,2};

% data_2_plot_2 = pw.stats_baseline{1,1};

steps         = diff(pw.full_trial.freq_baseline); % according to the fft time window
freq2plot_1   = 2:steps(1):12;
closestfreq_1 = dsearchn(pw.full_trial.freq_baseline,freq2plot_1');
freq_vector   = pw.full_trial.freq_baseline(closestfreq_1);

% Choose channels to plot
% Totty data -> 1-7   mPFC PL
%            -> 8-14  mPFC IL
%            -> 15-29 HPC

% Monopolar - common channels for all animals after preprocessing
channels{1} = 1:8;
channels{2} = 9:1;
channels{3} = 14:21;


% % Bipolar
% channels{1} = 1:5;
% channels{2} = 6:10;
% channels{3} = 11:19;

% Number of CS-trials
%CSIT = (1:9)'*(1:5); % Extinction
%CSIT = (1:5);
CSIT = 1;

figure('WindowState','maximized');
set(gcf,'color','w');

for ii = 1:length(channels)
    for jj = 1:size(CSIT,1)


        sgtitle({'Welch power spectral density estimate ';['(window = ' num2str(pw.full_trial.baseline_timewin./1000) 's' ' - ' 'overlap = ' num2str(pw.full_trial.overlap) '%)']})
        set(gcf,'color','white')

        if ii == 2
            pp = jj + size(CSIT,1);
        elseif ii == 3
            pp = jj + size(CSIT,1)*2;
        else
            pp = jj;
        end

        subplot(3,9,pp)
        plot(freq_vector, mean(data_2_plot_1(channels{ii},:),1) ,'linew', 2,'Color',[.6, .6, .6]);
        hold on
        plot(freq_vector, mean(mean(data_2_plot_2(channels{ii},:,CSIT(jj)),3),1) ,'linew', 2,'Color',[.6, 0, 0]);
        xline(freq_vector(132))

        xlabel('(Hz)','FontSize',9), ylabel('Normalized Power (A.U.)','FontSize',11)
        xlim([1 12])

        if ii == 1
            ylim([0 0.015])
        elseif ii == 2
            ylim([0 0.015])
        else
            ylim([0 0.015])
        end

        title(['CS Blk: ' num2str(jj)])



    end
end


% Settings
% newStr1 = id(1:end-20);
% name = strcat(Path,'/',newStr1,'_pw_total_power_norm');

% % save figure
%saveas(gcf,name,'png')

% close all

clear('newStr1','channels','steps','CSIT','closestfreq_1','closestfreq','freq_vector','data_2_plot_1','data_2_plot_2','freq2plot_1','ii','jj','name','path','pp') 


%% Plot to check 5 block session average. Averaging Channels per substrate as well
% Baseline normalization

%data_2_plot_1 = pw.stats_total_power{1,1};
data_2_plot_2 = pw.stats_baseline{1,1};

% data_2_plot_2 = pw.stats_baseline{1,1};

steps         = diff(pw.full_trial.freq_baseline); % according to the fft time window
freq2plot_1   = 2:steps(1):12;
closestfreq_1 = dsearchn(pw.full_trial.freq_baseline,freq2plot_1');
freq_vector   = pw.full_trial.freq_baseline(closestfreq_1);

% Choose channels to plot
% Totty data -> 1-7   mPFC PL
%            -> 8-14  mPFC IL
%            -> 15-29 HPC

% Monopolar - common channels for all animals after preprocessing
channels{1} = 1;
channels{2} = 7:11;
channels{3} = 14:21;


% % Bipolar
% channels{1} = 1:5;
% channels{2} = 6:10;
% channels{3} = 11:19;

% Number of CS-trials
CSIT = (1:9)'*(1:5); % Extinction
%CSIT = (1:5);


figure('WindowState','maximized');
set(gcf,'color','w');

for ii = 1:length(channels)
    for jj = 1:size(CSIT,1)


        sgtitle({'Welch power spectral density estimate ';['(window = ' num2str(pw.full_trial.baseline_timewin./1000) 's' ' - ' 'overlap = ' num2str(pw.full_trial.overlap) '%)']})
        set(gcf,'color','white')

        if ii == 2
            pp = jj + size(CSIT,1);
        elseif ii == 3
            pp = jj + size(CSIT,1)*2;
        else
            pp = jj;
        end

        subplot(3,9,pp)
%         plot(freq_vector, mean(data_2_plot_1(channels{ii},:),1) ,'linew', 2,'Color',[.6, .6, .6]);
%         hold on
        plot(freq_vector, mean(mean(data_2_plot_2(channels{ii},:,CSIT(jj)),3),1) ,'linew', 2,'Color',[.6, 0, 0]);
        xline(freq_vector(132))

        xlabel('(Hz)','FontSize',9), ylabel('Normalized Power (A.U.)','FontSize',11)
        xlim([1 12])

        if ii == 1
            ylim([0 2.5])
        elseif ii == 2
            ylim([0 2.5])
        else
            ylim([0 2.5])
        end
        

        title(['CS Blk: ' num2str(jj)])



    end
end

clear('','','','')


% % Settings
% newStr1 = id(1:end-20);
% name = strcat(Path,'/',newStr1,'_pw_baseline_norm');
% 
% % save figure
% saveas(gcf,name,'png')
% 
% close all

clear('newStr1','channels','steps','CSIT','closestfreq_1','closestfreq','freq_vector','data_2_plot_1','data_2_plot_2','freq2plot_1','ii','jj','name','path','pp') 

%% Plot to check full session. Channels per substrate or Average channels. ---> CHECK THE CODE !!! function mean

data_2_plot_1 = pw.stats_total_power{1,1};
%data_2_plot_2 = pw.stats_baseline{1,1};
data_2_plot_2 = pw.stats_total_power{1,2};

steps         = diff(pw.full_trial.freq_baseline); % according to the fft time window
freq2plot_1   = 2:steps(1):12;
closestfreq_1 = dsearchn(pw.full_trial.freq_baseline,freq2plot_1');
freq_vector   = pw.full_trial.freq_baseline(closestfreq_1);

% Choose channels to plot
% Tugce data -> 1-16   RE

channels = size(data_2_plot_2,1);

% Number of CS-trials and ITI to average
CSIT = 1:4

figure('WindowState','maximized');

% for ii = 1:length(channels)


    %subplot(4,4,ii)

    sgtitle({'Welch power spectral density estimate ';['(window = ' num2str(pw.full_trial.baseline_timewin./1000) 's' ' - ' 'overlap = ' num2str(pw.full_trial.overlap) '%)']}) 
    set(gcf,'color','white')

    %plot(freq_vector, mean(data_2_plot_1(channels(ii),:),3) ,'linew', 2,'Color',[.6, .6, .6]);
    %plot(freq_vector, mean(mean(data_2_plot_1(channels,:),3),1) ,'linew', 2,'Color',[.6, .6, .6]);
    hold on
    
    for jj = 1:length(CSIT)
        %plot(freq_vector, mean(data_2_plot_2(channels(ii),:),3) ,'linew', 2,'Color',[.6, .6, .6]);
        plot(freq_vector, mean(mean(data_2_plot_2(channels,:,CSIT(jj)),3),1) ,'linew', 2);
    end

    xlabel('(Hz)','FontSize',9), ylabel('Normalized Power (A.U.)','FontSize',11)
    %xlim([freq2plot(1) freq2plot(end)])
    xlim([3 10])
    %ylim([0 0.012])
    %ylim([0 2.5])

    %title(['Channel: ' num2str(channels(ii))])

%     if ii == length(channels)
        %legend ('baseline -> 52%','CS-1 -> 44%','CS-2 -> 92%','CS-3 -> 34%','CS-4 -> 26%','CS-5 -> 0%')
        legend ('CS-1 -> 44%','CS-2 -> 92%','CS-3 -> 34%','CS-4 -> 26%','CS-5 -> 0%')

%     end
  

% end


clear('ch','closestfreq','CSIT','freq2plot','ii','steps', 'total_p_baseline','total_p_CS_trial','total_p_ITI')

%% FOOOF - fitting oscillations & one over f
% %          Spectral Parameterization - Matlab Wrapper
% 
% %          >> https://github.com/fooof-tools/fooof
% 
% 
% % First: Create an environment in Anaconda for python 3.8.
% %        (Latest update for FOOOF uses functions for version 3.8)
% %        Check what is the Matlab limitation. 2022a set up to Python 3.8.
% 
% 
% % FOOOF settings
% settings = struct();  % Defaults values
% pw.over_f.f_range = [1, 50];
% return_model = true;
% 
% % Baseline
% % Preallocate - one channel per row
% pw.over_f.baseline = struct('Analise', cell(1, size(pw.full_trial.Pxx{1,1},1)),'NormOriginal_fit', cell(1, size(pw.full_trial.Pxx{1,1},1)),'NormModel_fit', cell(1, size(pw.full_trial.Pxx{1,1},1)));
% 
% for ff = 1:size(pw.full_trial.Pxx{1,1},1)
%     
%     pw.over_f.baseline(ff).Analise = fooof(pw.full_trial.freq_baseline, pw.full_trial.Pxx{1,1}(ff,:)', pw.over_f.f_range, settings, return_model); % Transpose PSD, to make inputs with the same direction
%        
%     % Normalize Spectrum
%     pw.over_f.baseline(ff).NormOriginal_fit = pw.over_f.baseline(ff).Analise.power_spectrum   - pw.over_f.baseline(ff).Analise.ap_fit; % Original Spectrum - Aperiodic Fit
%     pw.over_f.baseline(ff).NormModel_fit    = pw.over_f.baseline(ff).Analise.fooofed_spectrum - pw.over_f.baseline(ff).Analise.ap_fit; % Full Model Fit - Aperiodic Fit
% end
% 
% 
% % CS-Trials
% % Preallocate - one channel per row
% pw.over_f.CS_Trials = struct('Analise', cell(size(pw.full_trial.Pxx{2,1},1), size(pw.full_trial.Pxx{2,1},3)),'NormOriginal_fit', cell(size(pw.full_trial.Pxx{2,1},1), size(pw.full_trial.Pxx{2,1},3)),'NormModel_fit', cell(size(pw.full_trial.Pxx{2,1},1), size(pw.full_trial.Pxx{2,1},3)));
% 
% for gg = 1:size(pw.full_trial.Pxx{2,1},3)
%     for ff = 1:size(pw.full_trial.Pxx{2,1},1)
% 
%         pw.over_f.CS_Trials(ff,gg).Analise = fooof(pw.full_trial.freq_CS_trials, pw.full_trial.Pxx{2,1}(ff,:,gg)', pw.over_f.f_range, settings, return_model); % Transpose PSD, to make inputs with the same direction
% 
%         % Normalize Spectrum
%         pw.over_f.CS_Trials(ff,gg).NormOriginal_fit = pw.over_f.CS_Trials(ff).Analise.power_spectrum   - pw.over_f.CS_Trials(ff).Analise.ap_fit; % Original Spectrum - Aperiodic Fit
%         pw.over_f.CS_Trials(ff,gg).NormModel_fit    = pw.over_f.CS_Trials(ff).Analise.fooofed_spectrum - pw.over_f.CS_Trials(ff).Analise.ap_fit; % Full Model Fit - Aperiodic Fit
%     end
% end
% 
% 
% % ITI
% % Preallocate - one channel per row
% pw.over_f.ITI = struct('Analise', cell(size(pw.full_trial.Pxx{3,1},1), size(pw.full_trial.Pxx{3,1},3)),'NormOriginal_fit', cell(size(pw.full_trial.Pxx{3,1},1), size(pw.full_trial.Pxx{3,1},3)),'NormModel_fit', cell(size(pw.full_trial.Pxx{3,1},1), size(pw.full_trial.Pxx{3,1},3)));
% 
% for gg = 1:size(pw.full_trial.Pxx{3,1},3)
%     for ff = 1:size(pw.full_trial.Pxx{3,1},1)
% 
%         pw.over_f.ITI(ff,gg).Analise = fooof(pw.full_trial.freq_ITI, pw.full_trial.Pxx{3,1}(ff,:,gg)', pw.over_f.f_range, settings, return_model); % Transpose PSD, to make inputs with the same direction
% 
%         % Normalize Spectrum
%         pw.over_f.ITI(ff,gg).NormOriginal_fit = pw.over_f.ITI(ff).Analise.power_spectrum   - pw.over_f.ITI(ff).Analise.ap_fit; % Original Spectrum - Aperiodic Fit
%         pw.over_f.ITI(ff,gg).NormModel_fit    = pw.over_f.ITI(ff).Analise.fooofed_spectrum - pw.over_f.ITI(ff).Analise.ap_fit; % Full Model Fit - Aperiodic Fit
%     end
% end
% 
% 
% % FOOF function to plot all analyse
% %fooof_plot(pw.over_f.baseline(1).Analise)
% 
% clear('pe','ff','return_model','settings')


%% Save

% Settings
newStr1 = id(1:end-20);
name = strcat(Path,'/',newStr1,'_pw');

% Save data
save(name,'pw','-v7.3')

clear('name','newStr1','path') 

%% last update 26/01/2024
%  listening: Sonic Youth - Disconnection Notice
