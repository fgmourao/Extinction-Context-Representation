% FOOOF - fitting oscillations & one over f
%          Spectral Parameterization - Matlab Wrapper

%          >> https://github.com/fooof-tools/fooof

% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  12/2023
% Last update: 02/2024

%%
% First: Create an environment in Anaconda for python 3.8.
%        (Latest update for FOOOF uses functions for version 3.8)
%        Check what is the Matlab limitation. 2022a set up to Python 3.8.

%pyenv(Version="~/anaconda3/envs/Python_3_10/bin/python")

% FOOOF settings
settings = struct();  % Defaults values
pw.over_f.f_range = [1, 40];
return_model = true;

% Baseline
% Preallocate - one channel per row
pw.over_f.baseline = struct('Analise', cell(1, size(pw.full_trial.Pxx{1,1},1)),'NormOriginal_fit', cell(1, size(pw.full_trial.Pxx{1,1},1)),'NormModel_fit', cell(1, size(pw.full_trial.Pxx{1,1},1)));

for ff = 1:size(pw.full_trial.Pxx{1,1},1)
    
    pw.over_f.baseline(ff).Analise = fooof(pw.full_trial.freq_baseline, pw.full_trial.Pxx{1,1}(ff,:)', pw.over_f.f_range, settings, return_model); % Transpose PSD, to make inputs with the same direction
       
    % Normalize Spectrum
    pw.over_f.baseline(ff).NormOriginal_fit = pw.over_f.baseline(ff).Analise.power_spectrum   - pw.over_f.baseline(ff).Analise.ap_fit; % Original Spectrum - Aperiodic Fit
    pw.over_f.baseline(ff).NormModel_fit    = pw.over_f.baseline(ff).Analise.fooofed_spectrum - pw.over_f.baseline(ff).Analise.ap_fit; % Full Model Fit - Aperiodic Fit
end


% CS-Trials
% Preallocate - one channel per row
pw.over_f.CS_Trials = struct('Analise', cell(size(pw.full_trial.Pxx{2,1},1), size(pw.full_trial.Pxx{2,1},3)),'NormOriginal_fit', cell(size(pw.full_trial.Pxx{2,1},1), size(pw.full_trial.Pxx{2,1},3)),'NormModel_fit', cell(size(pw.full_trial.Pxx{2,1},1), size(pw.full_trial.Pxx{2,1},3)));


for gg = 1:length(CSIT{ms}) % gg = 1:size(pw.full_trial.Pxx{2,1},3)
    for ff = 1:size(pw.full_trial.Pxx{2,1},1)

        pw.over_f.CS_Trials(ff,gg).Analise = fooof(pw.full_trial.freq_CS_trials, pw.full_trial.Pxx{2,1}(ff,:,gg)', pw.over_f.f_range, settings, return_model); % Transpose PSD, to make inputs with the same direction

        % Normalize Spectrum
        pw.over_f.CS_Trials(ff,gg).NormOriginal_fit = pw.over_f.CS_Trials(ff,gg).Analise.power_spectrum   - pw.over_f.CS_Trials(ff,gg).Analise.ap_fit; % Original Spectrum - Aperiodic Fit
        pw.over_f.CS_Trials(ff,gg).NormModel_fit    = pw.over_f.CS_Trials(ff,gg).Analise.fooofed_spectrum - pw.over_f.CS_Trials(ff,gg).Analise.ap_fit; % Full Model Fit - Aperiodic Fit
    end
end


% ITI
% Preallocate - one channel per row
pw.over_f.ITI = struct('Analise', cell(size(pw.full_trial.Pxx{3,1},1), size(pw.full_trial.Pxx{3,1},3)),'NormOriginal_fit', cell(size(pw.full_trial.Pxx{3,1},1), size(pw.full_trial.Pxx{3,1},3)),'NormModel_fit', cell(size(pw.full_trial.Pxx{3,1},1), size(pw.full_trial.Pxx{3,1},3)));

for gg = 1:length(CSIT{ms}) %gg = 1:size(pw.full_trial.Pxx{3,1},3)
    for ff = 1:size(pw.full_trial.Pxx{3,1},1)

        pw.over_f.ITI(ff,gg).Analise = fooof(pw.full_trial.freq_ITI, pw.full_trial.Pxx{3,1}(ff,:,gg)', pw.over_f.f_range, settings, return_model); % Transpose PSD, to make inputs with the same direction

        % Normalize Spectrum
        pw.over_f.ITI(ff,gg).NormOriginal_fit = pw.over_f.ITI(ff,gg).Analise.power_spectrum   - pw.over_f.ITI(ff).Analise.ap_fit; % Original Spectrum - Aperiodic Fit
        pw.over_f.ITI(ff,gg).NormModel_fit    = pw.over_f.ITI(ff,gg).Analise.fooofed_spectrum - pw.over_f.ITI(ff).Analise.ap_fit; % Full Model Fit - Aperiodic Fit
    end
end


% FOOF function to plot all analyse
% fooof_plot(pw.over_f.CS_Trials(1, 3).Analise)

clear('pe','ff','return_model','settings')

%% Extract data

theta_range = [2 12];
closestfreq_1  = find(theta_range(1)<pw.over_f.CS_Trials(1, 1).Analise.freqs & pw.over_f.CS_Trials(1, 1).Analise.freqs<theta_range(2));
freq1_1 = pw.over_f.CS_Trials(1, 1).Analise.freqs(1,closestfreq_1);

% full spectrum

for jj = 1:size(pw.over_f.CS_Trials,1)
    for ii = 1:length(CSIT{ms})

        pw.over_f.stats_total_power{1,1}(jj,:) = pw.over_f.baseline(jj).NormModel_fit(1,closestfreq_1);   % Baseline
        pw.over_f.stats_total_power{1,2}(jj,:,ii) = pw.over_f.CS_Trials(jj, ii).NormModel_fit(1,closestfreq_1); % CS-Trials
        pw.over_f.stats_total_power{1,3}(jj,:,ii) = pw.over_f.ITI(jj, ii).NormModel_fit(1,closestfreq_1); % ITI

    end
end

% Averaging Trials
pw.over_f.stats_total_power_mean{1,1} = 10.*(pw.over_f.stats_total_power{1,1}./sum(squeeze(mean(pw.over_f.stats_total_power{1,1},3)),2));                   % Baseline
pw.over_f.stats_total_power_mean{1,2} = 10.*(squeeze(mean(pw.over_f.stats_total_power{1,2},3))./sum(squeeze(mean(pw.over_f.stats_total_power{1,2},3)),2));   % CS-Trials
pw.over_f.stats_total_power_mean{1,3} = 10.*(squeeze(mean(pw.over_f.stats_total_power{1,3},3))./sum(squeeze(mean(pw.over_f.stats_total_power{1,2},3)),2));   % ITI



range = [3, 6];
closestfreq_2  = find(range(1)<freq1_1 & freq1_1<range(2));
freq2_2 = freq1_1 (1,closestfreq_2);

% Averaging band frequency
pw.over_f.stats_total_power_mean{2,1} = mean(pw.over_f.stats_total_power_mean{1,1}(:,closestfreq_2),2);          % Baseline
pw.over_f.stats_total_power_mean{2,2} = squeeze(mean(pw.over_f.stats_total_power_mean{1,2}(:,closestfreq_2),2)); % CS-Trials
pw.over_f.stats_total_power_mean{2,3} = squeeze(mean(pw.over_f.stats_total_power_mean{1,3}(:,closestfreq_2),2)); % ITI

range = [3.6, 4.2];
closestfreq_2  = find(range(1)<freq1_1 & freq1_1<range(2));
freq2_2 = freq1_1 (1,closestfreq_2);

% Peak frequency
pw.over_f.stats_total_power_peak{2,1} = max(pw.over_f.stats_total_power_mean{1,1}(:,closestfreq_2),[],2);          % Baseline
pw.over_f.stats_total_power_peak{2,2} = squeeze(max(pw.over_f.stats_total_power_mean{1,2}(:,closestfreq_2),[],2)); % CS-Trials
pw.over_f.stats_total_power_peak{2,3} = squeeze(max(pw.over_f.stats_total_power_mean{1,3}(:,closestfreq_2),[],2)); % ITI


% 6 - 9 Hz
range = [6, 9];
closestfreq_3  = find(range(1)<freq1_1 & freq1_1<range(2));
freq3_3 = freq1_1 (1,closestfreq_3);

% Averaging band frequency
pw.over_f.stats_total_power_mean{3,1} = mean(pw.over_f.stats_total_power_mean{1,1}(:,closestfreq_3),2);          % Baseline
pw.over_f.stats_total_power_mean{3,2} = squeeze(mean(pw.over_f.stats_total_power_mean{1,2}(:,closestfreq_3),2)); % CS-Trials
pw.over_f.stats_total_power_mean{3,3} = squeeze(mean(pw.over_f.stats_total_power_mean{1,3}(:,closestfreq_3),2)); % ITI

range = [6.5, 7.2];
closestfreq_3  = find(range(1)<freq1_1 & freq1_1<range(2));
freq2_2 = freq1_1 (1,closestfreq_2);

% peak frequency
pw.over_f.stats_total_power_peak{3,1} = max(pw.over_f.stats_total_power_mean{1,1}(:,closestfreq_3),[],2);          % Baseline
pw.over_f.stats_total_power_peak{3,2} = squeeze(max(pw.over_f.stats_total_power_mean{1,2}(:,closestfreq_3),[],2)); % CS-Trials
pw.over_f.stats_total_power_peak{3,3} = squeeze(max(pw.over_f.stats_total_power_mean{1,3}(:,closestfreq_3),[],2)); % ITI

%% Save

% Settings
%ms = 1;
Path    = files.FilesLoaded{1,1}(ms).folder;
newStr1 = id(1:end-8);
name_1 = strcat(Path,'/',newStr1,'_pw_foof_nFFT_32768');


% Save data
save(name_1,'pw','-v7.3')

clear('name','newStr1','path') 

%%
