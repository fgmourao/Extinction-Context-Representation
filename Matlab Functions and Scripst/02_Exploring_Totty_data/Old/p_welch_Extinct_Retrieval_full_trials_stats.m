
%%  Welch power spectral density estimate - Stats

% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  12/2023
% Last update: 06/2024


%% Extract data from theta range (2 - 12 Hz). Total power normalization

% Set freuency6 range
pw.full_trial.parameters.frex_3_6Hz     = 3:steps(1):6;
pw.full_trial.parameters.frex_idx_3_6Hz = dsearchn(pw.full_trial.parameters.frex_2_12Hz',pw.full_trial.parameters.frex_3_6Hz');

pw.full_trial.parameters.frex_6_8Hz     = 6:steps(1):8;
pw.full_trial.parameters.frex_idx_6_8Hz = dsearchn(pw.full_trial.parameters.frex_2_12Hz',pw.full_trial.parameters.frex_6_8Hz');

pw.full_trial.parameters.frex_7_10Hz     = 7:steps(1):10;
pw.full_trial.parameters.frex_idx_7_10Hz = dsearchn(pw.full_trial.parameters.frex_2_12Hz',pw.full_trial.parameters.frex_7_10Hz');


% Integrated band frequencies
% Integrated band frequency - All Trials
pw.full_trial.Pxx_TotalPower_norm{1,2} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm{1,1}(:,pw.full_trial.parameters.frex_idx_3_6Hz,:),2)); % Baseline
pw.full_trial.Pxx_TotalPower_norm{2,2} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm{2,1}(:,pw.full_trial.parameters.frex_idx_3_6Hz,:),2)); % CS-Trials
pw.full_trial.Pxx_TotalPower_norm{3,2} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm{3,1}(:,pw.full_trial.parameters.frex_idx_3_6Hz,:),2)); % ITI

pw.full_trial.Pxx_TotalPower_norm{1,3} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm{1,1}(:,pw.full_trial.parameters.frex_idx_6_8Hz,:),2)); % Baseline
pw.full_trial.Pxx_TotalPower_norm{2,3} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm{2,1}(:,pw.full_trial.parameters.frex_idx_6_8Hz,:),2)); % CS-Trials
pw.full_trial.Pxx_TotalPower_norm{3,3} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm{3,1}(:,pw.full_trial.parameters.frex_idx_6_8Hz,:),2)); % ITI

pw.full_trial.Pxx_TotalPower_norm{1,4} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm{1,1}(:,pw.full_trial.parameters.frex_idx_7_10Hz,:),2)); % Baseline
pw.full_trial.Pxx_TotalPower_norm{2,4} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm{2,1}(:,pw.full_trial.parameters.frex_idx_7_10Hz,:),2)); % CS-Trials
pw.full_trial.Pxx_TotalPower_norm{3,4} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm{3,1}(:,pw.full_trial.parameters.frex_idx_7_10Hz,:),2)); % ITI


% Integrated band frequency - First 10 Trials
pw.full_trial.Pxx_TotalPower_norm_mean_first_trials{1,2} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm_mean_first_trials{1,1}(:,pw.full_trial.parameters.frex_idx_3_6Hz,:),2)); % Baseline
pw.full_trial.Pxx_TotalPower_norm_mean_first_trials{2,2} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm_mean_first_trials{2,1}(:,pw.full_trial.parameters.frex_idx_3_6Hz,:),2)); % CS-Trials
pw.full_trial.Pxx_TotalPower_norm_mean_first_trials{3,2} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm_mean_first_trials{3,1}(:,pw.full_trial.parameters.frex_idx_3_6Hz,:),2)); % ITI

pw.full_trial.Pxx_TotalPower_norm_mean_first_trials{1,3} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm_mean_first_trials{1,1}(:,pw.full_trial.parameters.frex_idx_6_8Hz,:),2)); % Baseline
pw.full_trial.Pxx_TotalPower_norm_mean_first_trials{2,3} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm_mean_first_trials{2,1}(:,pw.full_trial.parameters.frex_idx_6_8Hz,:),2)); % CS-Trials
pw.full_trial.Pxx_TotalPower_norm_mean_first_trials{3,3} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm_mean_first_trials{3,1}(:,pw.full_trial.parameters.frex_idx_6_8Hz,:),2)); % ITI

pw.full_trial.Pxx_TotalPower_norm_mean_first_trials{1,4} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm_mean_first_trials{1,1}(:,pw.full_trial.parameters.frex_idx_7_10Hz,:),2)); % Baseline
pw.full_trial.Pxx_TotalPower_norm_mean_first_trials{2,4} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm_mean_first_trials{2,1}(:,pw.full_trial.parameters.frex_idx_7_10Hz,:),2)); % CS-Trials
pw.full_trial.Pxx_TotalPower_norm_mean_first_trials{3,4} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm_mean_first_trials{3,1}(:,pw.full_trial.parameters.frex_idx_7_10Hz,:),2)); % ITI


% Integrated band frequency - Last 10 Trials
pw.full_trial.Pxx_TotalPower_norm_mean_last_trials{1,2} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm_mean_last_trials{1,1}(:,pw.full_trial.parameters.frex_idx_3_6Hz,:),2)); % Baseline
pw.full_trial.Pxx_TotalPower_norm_mean_last_trials{2,2} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm_mean_last_trials{2,1}(:,pw.full_trial.parameters.frex_idx_3_6Hz,:),2)); % CS-Trials
pw.full_trial.Pxx_TotalPower_norm_mean_last_trials{3,2} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm_mean_last_trials{3,1}(:,pw.full_trial.parameters.frex_idx_3_6Hz,:),2)); % ITI

pw.full_trial.Pxx_TotalPower_norm_mean_last_trials{1,3} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm_mean_last_trials{1,1}(:,pw.full_trial.parameters.frex_idx_6_8Hz,:),2)); % Baseline
pw.full_trial.Pxx_TotalPower_norm_mean_last_trials{2,3} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm_mean_last_trials{2,1}(:,pw.full_trial.parameters.frex_idx_6_8Hz,:),2)); % CS-Trials
pw.full_trial.Pxx_TotalPower_norm_mean_last_trials{3,3} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm_mean_last_trials{3,1}(:,pw.full_trial.parameters.frex_idx_6_8Hz,:),2)); % ITI

pw.full_trial.Pxx_TotalPower_norm_mean_last_trials{1,4} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm_mean_last_trials{1,1}(:,pw.full_trial.parameters.frex_idx_7_10Hz,:),2)); % Baseline
pw.full_trial.Pxx_TotalPower_norm_mean_last_trials{2,4} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm_mean_last_trials{2,1}(:,pw.full_trial.parameters.frex_idx_7_10Hz,:),2)); % CS-Trials
pw.full_trial.Pxx_TotalPower_norm_mean_last_trials{3,4} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm_mean_last_trials{3,1}(:,pw.full_trial.parameters.frex_idx_7_10Hz,:),2)); % ITI



% Peak Frequency
% Peak frequency range - All Trials
pw.full_trial.Pxx_TotalPower_norm_peak{1,2} = squeeze(max(pw.full_trial.Pxx_TotalPower_norm{1,1}(:,pw.full_trial.parameters.frex_idx_3_6Hz,:),[],2)); % Baseline
pw.full_trial.Pxx_TotalPower_norm_peak{2,2} = squeeze(max(pw.full_trial.Pxx_TotalPower_norm{2,1}(:,pw.full_trial.parameters.frex_idx_3_6Hz,:),[],2)); % CS-Trials
pw.full_trial.Pxx_TotalPower_norm_peak{3,2} = squeeze(max(pw.full_trial.Pxx_TotalPower_norm{3,1}(:,pw.full_trial.parameters.frex_idx_3_6Hz,:),[],2)); % ITI

pw.full_trial.Pxx_TotalPower_norm_peak{1,3} = squeeze(max(pw.full_trial.Pxx_TotalPower_norm{1,1}(:,pw.full_trial.parameters.frex_idx_6_8Hz,:),[],2)); % Baseline
pw.full_trial.Pxx_TotalPower_norm_peak{2,3} = squeeze(max(pw.full_trial.Pxx_TotalPower_norm{2,1}(:,pw.full_trial.parameters.frex_idx_6_8Hz,:),[],2)); % CS-Trials
pw.full_trial.Pxx_TotalPower_norm_peak{3,3} = squeeze(max(pw.full_trial.Pxx_TotalPower_norm{3,1}(:,pw.full_trial.parameters.frex_idx_6_8Hz,:),[],2)); % ITI

pw.full_trial.Pxx_TotalPower_norm_peak{1,4} = squeeze(max(pw.full_trial.Pxx_TotalPower_norm{1,1}(:,pw.full_trial.parameters.frex_idx_7_10Hz,:),[],2)); % Baseline
pw.full_trial.Pxx_TotalPower_norm_peak{2,4} = squeeze(max(pw.full_trial.Pxx_TotalPower_norm{2,1}(:,pw.full_trial.parameters.frex_idx_7_10Hz,:),[],2)); % CS-Trials
pw.full_trial.Pxx_TotalPower_norm_peak{3,4} = squeeze(max(pw.full_trial.Pxx_TotalPower_norm{3,1}(:,pw.full_trial.parameters.frex_idx_7_10Hz,:),[],2)); % ITI


%% Save

%newStr = regexprep(files.id.name,'.mat','_');
%newStr = files.id(ms).name(1:end-8);
newStr = id(1:end-8);

path = '/Users/flavio/Desktop';
%path = files.FilesLoaded{1, 1}.folder;

name = strcat(path,'/',newStr,'_pw_Total_Power_win_1000_nFFT_32768_fullTrial_2048TimeW');

% save data
save(name,'pw','-v7.3')

clear('name','newStr','path')

%%
