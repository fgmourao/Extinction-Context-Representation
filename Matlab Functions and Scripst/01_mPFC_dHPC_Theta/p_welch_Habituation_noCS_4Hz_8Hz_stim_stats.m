
%%  Welch power spectral density estimate - Stats

% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  12/2023
% Last update: 01/2024


%% Extract data
% Baseline Normalization 

% Define band frequency
% 2 - 12 Hertz
steps         = diff(pw.full_trial.freq_baseline); % according to the fft time window
freq2plot     = 2:steps(1):12;
pw.stats_baseline{5,1}  = dsearchn(pw.full_trial.freq_baseline,freq2plot');

if ~isempty(pw.full_trial.Pxx{2,1})
    % CS-Trials
    pw.stats_baseline{1,1} = pw.full_trial.Pxx{2,1}(:,pw.stats_baseline{5,1},:)./pw.full_trial.Pxx{1,1}(:,pw.stats_baseline{5,1}); % / Baseline
    pw.stats_baseline{2,1} = 10*log10(pw.full_trial.Pxx{2,1}(:,pw.stats_baseline{5,1},:)./pw.full_trial.Pxx{1,1}(:,pw.stats_baseline{5,1})); % log power / Baseline
    pw.stats_baseline{3,1} = (pw.full_trial.Pxx{2,1}(:,pw.stats_baseline{5,1},:) - pw.full_trial.Pxx{1,1}(:,pw.stats_baseline{5,1}))./pw.full_trial.Pxx{1,1}(:,pw.stats_baseline{5,1}); % / percentage
    pw.stats_baseline{4,1} = (pw.full_trial.Pxx{2,1}(:,pw.stats_baseline{5,1},:) - mean(pw.full_trial.Pxx{1,1}(:,pw.stats_baseline{5,1}),2))./ std(pw.full_trial.Pxx{1,1}(:,pw.stats_baseline{5,1}),[],2); % / Z-score
end

clear('steps','freq2plot')

%%
% Total Power Normalization according to the Maren & Totty paper

% Define band frequency
% 2 - 12 Hertz / Keep fullspec to plot
% steps         = diff(pw.full_trial.freq_CS_trials); % according to the fft time window
% freq2plot_1   = 2:steps(1):12;
% closestfreq_1 = dsearchn(pw.full_trial.freq_CS_trials,freq2plot_1');

theta_range = [2 12];
closestfreq_1  = find(theta_range(1)<pw.full_trial.freq_baseline & pw.full_trial.freq_baseline<theta_range(2));
freq1_1 = pw.full_trial.freq_baseline(closestfreq_1,1);

range = [3, 6];
closestfreq_2  = find(range(1)<freq1_1 & freq1_1<range(2));
freq2_2 = freq1_1 (closestfreq_2,1);

range = [6, 9];
closestfreq_3  = find(range(1)<freq1_1 & freq1_1<range(2));
freq3_3 = freq1_1 (closestfreq_3,1);


if ~isempty(pw.full_trial.Pxx{1,1})
    pw.stats.total_power{1,1} = pw.full_trial.Pxx{1,1}(:,closestfreq_1)./sum(pw.full_trial.Pxx{1,1}(:,closestfreq_1),2); % Baseline
    pw.stats.total_power{2,1} = mean(pw.stats.total_power{1,1}(:,closestfreq_2),2); % Baseline
    pw.stats.total_power{3,1} = mean(pw.stats.total_power{1,1}(:,closestfreq_3),2); % Baseline

end

if ~isempty(pw.full_trial.Pxx{2,1})
    pw.stats.total_power{1,2} = pw.full_trial.Pxx{2,1}(:,closestfreq_1,:)./sum(pw.full_trial.Pxx{2,1}(:,closestfreq_1),2); % CS-Trials
    pw.stats.total_power{2,2} = squeeze(mean(pw.stats.total_power{1,2}(:,closestfreq_2,:),2)); % CS-Trials
    pw.stats.total_power{3,2} = squeeze(mean(pw.stats.total_power{1,2}(:,closestfreq_3,:),2)); % CS-Trials

end

if ~isempty(pw.full_trial.Pxx{3,1})
    pw.stats.total_power{1,3} = pw.full_trial.Pxx{3,1}(:,closestfreq_1,:)./sum(pw.full_trial.Pxx{3,1}(:,closestfreq_1),2); % ITI
    pw.stats.total_power{2,3} = squeeze(mean(pw.stats.total_power{1,3}(:,closestfreq_2,:),2)); % ITI
    pw.stats.total_power{3,3} = squeeze(mean(pw.stats.total_power{1,3}(:,closestfreq_3,:),2)); % ITI

end


% Peak
pw.stats.total_power_peak{2,1} = max(pw.stats.total_power{1,1}(:,closestfreq_2),[],2); % Baseline
pw.stats.total_power_peak{2,2} = squeeze(max(pw.stats.total_power{1,2}(:,closestfreq_2,:),[],2)); % CS-Trials
pw.stats.total_power_peak{2,3} = squeeze(max(pw.stats.total_power{1,3}(:,closestfreq_2,:),[],2)); % ITI

% Peak
pw.stats.total_power_peak{3,1} = max(pw.stats.total_power{1,1}(:,closestfreq_3),[],2); % Baseline
pw.stats.total_power_peak{3,2} = squeeze(max(pw.stats.total_power{1,2}(:,closestfreq_3,:),[],2)); % CS-Trials
pw.stats.total_power_peak{3,3} = squeeze(max(pw.stats.total_power{1,3}(:,closestfreq_3,:),[],2)); % ITI



% Following the laboratory method, averages were calculated for every 5 trials.
% Graph based on Maren`s paper - nature communication

pw.stats.total_power_mean{2,2} = mean(pw.stats.total_power{2,2}(:,CSIT{ms}),2);  % CS-Trials
pw.stats.total_power_mean{2,3} = mean(pw.stats.total_power{2,3}(:,CSIT_1{ms}),2); % ITI

pw.stats.total_power_mean{3,2} = mean(pw.stats.total_power{3,2}(:,CSIT{ms}),2);   % CS-Trials
pw.stats.total_power_mean{3,3} = mean(pw.stats.total_power{3,3}(:,CSIT_1{ms}),2); % ITI


% Peak
pw.stats.total_power_peak_mean{2,2} = mean(pw.stats.total_power_peak{2,2}(:,CSIT{ms}),2);
pw.stats.total_power_peak_mean{2,3} = mean(pw.stats.total_power_peak{2,3}(:,CSIT_1{ms}),2);

pw.stats.total_power_peak_mean{3,2} = mean(pw.stats.total_power_peak{3,2}(:,CSIT{ms}),2);
pw.stats.total_power_peak_mean{3,3} = mean(pw.stats.total_power_peak{3,3}(:,CSIT_1{ms}),2);


clear('steps','freq2plot_1','freq2plot_2','closestfreq_2','freq2plot_3','closestfreq_3')

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
