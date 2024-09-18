%% Organize Granger

% Plot spectral causal graph. Version 2.
% Plot pairwise quantities in |P|, a 3-dim numerical matrix with
% first index representing target ("to"), second index source ("from")
% and third index frequency range - typically spectral causalities

% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  02/2024
% Last update: 02/2024

% Organizing the matrices into a single variable. Each row represents the mean of the chosen trials in a respective pair of channels."
% Respective plots. Spectra and graphs

%% Pair channels combination
% Row 1: PL  --> IL
% Row 2: IL  --> PL
% Row 3: PL  --> HPC
% Row 4: HPC --> PL
% Row 5: IL  --> HPC
% Row 6: HPC --> IL

% All Possibles  combinations
combinations_ = zeros(size(mvgc.Fint_3_6Hz{1, 1},1)*2,2);
combinations_(1:2:6,:) = flip(nchoosek(1:size(mvgc.data,1),2),2);
combinations_(2:2:6,:) = nchoosek(1:size(mvgc.data,1),2);

%% Data to Stats on Prism

% Matrices Organized into a single variable. Each row represents the mean of the chosen trials in a respective pair of channels."
% Respective plots. Spectra and graphs

% Baseline, CS-Trials and ITI
% This for-loop needs improvement due to redundancy in the baseline extraction data

% Baseline
for jj = 1:size(combinations_,1)

    mvgc.stats.Fint_3_6Hz_baseline_2_plot(jj,1) = mvgc.Fint_3_6Hz{1,1}(combinations_(jj,1), combinations_(jj,2)); % baseline
    mvgc.stats.Fint_6_8Hz_baseline_2_plot(jj,1) = mvgc.Fint_6_8Hz{1,1}(combinations_(jj,1), combinations_(jj,2)); % baseline
    mvgc.stats.Fint_7_10Hz_baseline_2_plot(jj,1) = mvgc.Fint_7_10Hz{1,1}(combinations_(jj,1), combinations_(jj,2)); % baseline

end

%CS-tone
for ii = 1:length(CSIT{ms})
    for jj = 1:size(combinations_,1)

        mvgc.stats.Fint_3_6Hz_CS_2_plot(jj,ii) = mvgc.Fint_3_6Hz{2,CSIT{ms}(ii)}(combinations_(jj,1), combinations_(jj,2)); % CS-Trials
        mvgc.stats.Fint_6_8Hz_CS_2_plot(jj,ii) = mvgc.Fint_6_8Hz{2,CSIT{ms}(ii)}(combinations_(jj,1), combinations_(jj,2)); % CS-Trials
        mvgc.stats.Fint_7_10Hz_CS_2_plot(jj,ii) = mvgc.Fint_7_10Hz{2,CSIT{ms}(ii)}(combinations_(jj,1), combinations_(jj,2)); % CS-Trials

    end
end

% ITI
for ii = 1:length(CSIT_1{ms})
    for jj = 1:size(combinations_,1)

        mvgc.stats.Fint_3_6Hz_ITI_2_plot(jj,ii) = mvgc.Fint_3_6Hz{3,CSIT_1{ms}(ii)}(combinations_(jj,1), combinations_(jj,2)); % ITI
        mvgc.stats.Fint_6_8Hz_ITI_2_plot(jj,ii) = mvgc.Fint_6_8Hz{3,CSIT_1{ms}(ii)}(combinations_(jj,1), combinations_(jj,2)); % ITI
        mvgc.stats.Fint_7_10Hz_ITI_2_plot(jj,ii) = mvgc.Fint_7_10Hz{3,CSIT_1{ms}(ii)}(combinations_(jj,1), combinations_(jj,2)); % ITI

    end
end

% Post-tones
for jj = 1:size(combinations_,1)

    mvgc.stats.Fint_3_6Hz_post_tone_2_plot(jj,1) = mvgc.Fint_3_6Hz{3,end}(combinations_(jj,1), combinations_(jj,2)); % Post-Tone
    mvgc.stats.Fint_6_8Hz_post_tone_2_plot(jj,1) = mvgc.Fint_6_8Hz{3,end}(combinations_(jj,1), combinations_(jj,2)); % Post-Tone
    mvgc.stats.Fint_7_10Hz_post_tone_2_plot(jj,1) = mvgc.Fint_7_10Hz{3,end}(combinations_(jj,1), combinations_(jj,2)); % Post-Tone

end


% Mean 10 first trials
mvgc.stats.Fint_3_6Hz_2_plot_mean_first(:,1) = mvgc.stats.Fint_3_6Hz_baseline_2_plot;              % baseline - Redundant but helped in the code below
mvgc.stats.Fint_3_6Hz_2_plot_mean_first(:,2) = mean(mvgc.stats.Fint_3_6Hz_CS_2_plot(:,1:(size(mvgc.stats.Fint_3_6Hz_CS_2_plot,2)./2)),2);    % CS
mvgc.stats.Fint_3_6Hz_2_plot_mean_first(:,3) = mean(mvgc.stats.Fint_3_6Hz_ITI_2_plot(:,1:(size(mvgc.stats.Fint_3_6Hz_ITI_2_plot,2)./2)),2);   % ITI

mvgc.stats.Fint_6_8Hz_2_plot_mean_first(:,1) = mvgc.stats.Fint_6_8Hz_baseline_2_plot;              % baseline - Redundant but helped in the code below
mvgc.stats.Fint_6_8Hz_2_plot_mean_first(:,2) = mean(mvgc.stats.Fint_6_8Hz_CS_2_plot(:,1:(size(mvgc.stats.Fint_6_8Hz_CS_2_plot,2)./2)),2);    % CS
mvgc.stats.Fint_6_8Hz_2_plot_mean_first(:,3) = mean(mvgc.stats.Fint_6_8Hz_ITI_2_plot(:,1:(size(mvgc.stats.Fint_6_8Hz_ITI_2_plot,2)./2)),2);   % ITI

mvgc.stats.Fint_7_10Hz_2_plot_mean_first(:,1) = mvgc.stats.Fint_7_10Hz_baseline_2_plot;              % baseline - Redundant but helped in the code below
mvgc.stats.Fint_7_10Hz_2_plot_mean_first(:,2) = mean(mvgc.stats.Fint_7_10Hz_CS_2_plot(:,1:(size(mvgc.stats.Fint_7_10Hz_CS_2_plot,2)./2)),2);    % CS
mvgc.stats.Fint_7_10Hz_2_plot_mean_first(:,3) = mean(mvgc.stats.Fint_7_10Hz_ITI_2_plot(:,1:(size(mvgc.stats.Fint_7_10Hz_ITI_2_plot,2)./2)),2);   % ITI


% Mean 10 last trials
mvgc.stats.Fint_3_6Hz_2_plot_mean_last(:,1) = mvgc.stats.Fint_3_6Hz_baseline_2_plot;               % baseline - Redundant but helped in the code below
mvgc.stats.Fint_3_6Hz_2_plot_mean_last(:,2) = mean(mvgc.stats.Fint_3_6Hz_CS_2_plot(:,(size(mvgc.stats.Fint_3_6Hz_CS_2_plot,2)./2+1:size(mvgc.stats.Fint_3_6Hz_CS_2_plot,2))),2);    % CS
mvgc.stats.Fint_3_6Hz_2_plot_mean_last(:,3) = mean(mvgc.stats.Fint_3_6Hz_ITI_2_plot(:,:),2);   % ITI

mvgc.stats.Fint_6_8Hz_2_plot_mean_last(:,1) = mvgc.stats.Fint_6_8Hz_baseline_2_plot;               % baseline - Redundant but helped in the code below
mvgc.stats.Fint_6_8Hz_2_plot_mean_last(:,2) = mean(mvgc.stats.Fint_6_8Hz_CS_2_plot(:,(size(mvgc.stats.Fint_3_6Hz_CS_2_plot,2)./2+1:size(mvgc.stats.Fint_6_8Hz_CS_2_plot,2))),2);    % CS
mvgc.stats.Fint_6_8Hz_2_plot_mean_last(:,3) = mean(mvgc.stats.Fint_6_8Hz_ITI_2_plot(:,:),2);   % ITI

mvgc.stats.Fint_7_10Hz_2_plot_mean_last(:,1) = mvgc.stats.Fint_7_10Hz_baseline_2_plot;               % baseline - Redundant but helped in the code below
mvgc.stats.Fint_7_10Hz_2_plot_mean_last(:,2) = mean(mvgc.stats.Fint_7_10Hz_CS_2_plot(:,(size(mvgc.stats.Fint_3_6Hz_CS_2_plot,2)./2+1:size(mvgc.stats.Fint_7_10Hz_CS_2_plot,2))),2);    % CS
mvgc.stats.Fint_7_10Hz_2_plot_mean_last(:,3) = mean(mvgc.stats.Fint_7_10Hz_ITI_2_plot(:,:),2);   % ITI

%% Specs Baseline, CS-Trials and ITI and Post-tone
% This for-loop needs improvement due to redundancy in the baseline extraction data

% Baseline
for jj = 1:size(combinations_,1)
    mvgc.stats.Spect_2_12Hz_baseline_2_plot(jj,:) = mvgc.F_spect{1,1}(combinations_(jj,1), combinations_(jj,2),:);% baseline
end


% CS
for ii = 1:length(CSIT{ms})
    for jj = 1:size(combinations_,1)
        mvgc.stats.Spect_2_12Hz_CS_2_plot(jj,ii,:) = mvgc.F_spect{2,CSIT{ms}(ii)}(combinations_(jj,1), combinations_(jj,2),:); % CS-Trials
    end
end

% ITI
for ii = 1:length(CSIT_1{ms})
    for jj = 1:size(combinations_,1)
        mvgc.stats.Spect_2_12Hz_ITI_2_plot(jj,ii,:) = mvgc.F_spect{3,CSIT_1{ms}(ii)}(combinations_(jj,1), combinations_(jj,2),:); % ITI
    end
end

% Post-Tone
for jj = 1:size(combinations_,1)
    mvgc.stats.Spect_2_12Hz_post_tone_2_plot(jj,:) = squeeze(mvgc.F_spect{3,end}(combinations_(jj,1), combinations_(jj,2),:)); % Post-Tone
end


%% Mean trials

% Mean first 10 trials - CS and ITI
mvgc.stats.Spect_2_12Hz_2_plot_mean_first = squeeze(mean(mvgc.stats.Spect_2_12Hz_CS_2_plot(:,1:(size(mvgc.stats.Fint_3_6Hz_CS_2_plot,2)./2),:),2));    % CS
%mvgc.stats.Spect_2_12Hz_2_plot_mean_first{2,1} = squeeze(mean(mvgc.stats.Spect_2_12Hz_ITI_2_plot(:,1:size(mvgc.stats.Fint_3_6Hz_CS_2_plot,2)./2,:),2));   % ITI


% Mean last 10 trials - CS and ITI
mvgc.stats.Spect_2_12Hz_2_plot_mean_last = squeeze(mean(mvgc.stats.Spect_2_12Hz_CS_2_plot(:,(size(mvgc.stats.Fint_3_6Hz_CS_2_plot,2)./2+1:size(mvgc.stats.Fint_3_6Hz_CS_2_plot,2)),:),2));    % CS
%mvgc.stats.Spect_2_12Hz_2_plot_mean_last(:,3,:) = squeeze(mean(mvgc.stats.Spect_2_12Hz_ITI_2_plot(:,1:size(mvgc.stats.Fint_3_6Hz_CS_2_plot,2)./2,:),2));   % ITI


%% peak frequency in range

% baseline 
mvgc.stats.Fint_3_6Hz_2_plot_peak_baseline  = max(mvgc.stats.Spect_2_12Hz_baseline_2_plot(:,mvgc.parameters.frex_idx_3_6Hz),[],2);
mvgc.stats.Fint_6_8Hz_2_plot_peak_baseline  = max(mvgc.stats.Spect_2_12Hz_baseline_2_plot(:,mvgc.parameters.frex_idx_6_8Hz),[],2);
mvgc.stats.Fint_7_10Hz_2_plot_peak_baseline = max(mvgc.stats.Spect_2_12Hz_baseline_2_plot(:,mvgc.parameters.frex_idx_7_10Hz),[],2);

% CS and ITI - peak over mean first 10 trials
mvgc.stats.Fint_3_6Hz_2_plot_peak_mean_first_trials  = max(mvgc.stats.Spect_2_12Hz_2_plot_mean_first(:,mvgc.parameters.frex_idx_3_6Hz),[],2);
mvgc.stats.Fint_6_8Hz_2_plot_peak_mean_first_trials  = max(mvgc.stats.Spect_2_12Hz_2_plot_mean_first(:,mvgc.parameters.frex_idx_6_8Hz),[],2);
mvgc.stats.Fint_7_10Hz_2_plot_peak_mean_first_trials = max(mvgc.stats.Spect_2_12Hz_2_plot_mean_first(:,mvgc.parameters.frex_idx_7_10Hz),[],2);

% CS and ITI - peak over mean last 10 trials
mvgc.stats.Fint_3_6Hz_2_plot_peak_mean_last_trials  = max(mvgc.stats.Spect_2_12Hz_2_plot_mean_last(:,mvgc.parameters.frex_idx_3_6Hz),[],3);
mvgc.stats.Fint_6_8Hz_2_plot_peak_mean_last_trials  = max(mvgc.stats.Spect_2_12Hz_2_plot_mean_last(:,mvgc.parameters.frex_idx_6_8Hz),[],3);
mvgc.stats.Fint_7_10Hz_2_plot_peak_mean_last_trials = max(mvgc.stats.Spect_2_12Hz_2_plot_mean_last(:,mvgc.parameters.frex_idx_7_10Hz),[],3);

% CS and ITI  -peak for each trial
mvgc.stats.Fint_3_6Hz_2_plot_peak_each_trials_CS   = max(mvgc.stats.Spect_2_12Hz_CS_2_plot(:,:,mvgc.parameters.frex_idx_3_6Hz),[],3);
mvgc.stats.Fint_6_8Hz_2_plot_peak_each_trials_CS   = max(mvgc.stats.Spect_2_12Hz_CS_2_plot(:,:,mvgc.parameters.frex_idx_6_8Hz),[],3);
mvgc.stats.Fint_3_6Hz_2_plot_peak_each_trials_ITI  = max(mvgc.stats.Spect_2_12Hz_ITI_2_plot(:,:,mvgc.parameters.frex_idx_3_6Hz),[],3);
mvgc.stats.Fint_6_8Hz_2_plot_peak_each_trials_ITI  = max(mvgc.stats.Spect_2_12Hz_ITI_2_plot(:,:,mvgc.parameters.frex_idx_6_8Hz),[],3);
mvgc.stats.Fint_7_10Hz_2_plot_peak_each_trials_ITI = max(mvgc.stats.Spect_2_12Hz_ITI_2_plot(:,:,mvgc.parameters.frex_idx_7_10Hz),[],3);
mvgc.stats.Fint_7_10Hz_2_plot_peak_each_trials_ITI = max(mvgc.stats.Spect_2_12Hz_ITI_2_plot(:,:,mvgc.parameters.frex_idx_7_10Hz),[],3);

% post tone 
mvgc.stats.Fint_3_6Hz_2_plot_peak_post_tone  = max(mvgc.stats.Spect_2_12Hz_post_tone_2_plot(:,mvgc.parameters.frex_idx_3_6Hz),[],2);
mvgc.stats.Fint_6_8Hz_2_plot_peak_post_tone  = max(mvgc.stats.Spect_2_12Hz_post_tone_2_plot(:,mvgc.parameters.frex_idx_6_8Hz),[],2);
mvgc.stats.Fint_7_10Hz_2_plot_peak_post_tone = max(mvgc.stats.Spect_2_12Hz_post_tone_2_plot(:,mvgc.parameters.frex_idx_7_10Hz),[],2);



%% Full width at half maximum
% for ii = 1:size(mvgc.stats.Spect_2_12Hz_2_plot_mean,1)
%     for jj = 1:size(mvgc.stats.Spect_2_12Hz_2_plot_mean,2)
% 
%         [mvgc.stats.FWHM_3_6Hz_samples(ii,jj), mvgc.stats.FWHM_3_6Hz_Hz(ii,jj)] = FWHM(squeeze(mvgc.stats.Spect_2_12Hz_2_plot_mean(ii,jj,mvgc.parameters.frex_idx_3_6Hz))',mvgc.parameters.frex_3_6Hz);
%         [mvgc.stats.FWHM_6_8Hz_samples(ii,jj), mvgc.stats.FWHM_6_8Hz_Hz(ii,jj)] = FWHM(squeeze(mvgc.stats.Spect_2_12Hz_2_plot_mean(ii,jj,mvgc.parameters.frex_idx_6_8Hz))',mvgc.parameters.frex_6_8Hz);
%    
%     end
% end

%% Save data

% Settings
%ms = 1;
Path    = files.FilesLoaded{1,1}(ms).folder;
newStr1 = id(1:end-8);
name_2 = strcat(Path,'/',newStr1,'_MVGC_Granger_with_data_to_prism');

% Save data
save(name_2,'mvgc','-v7.3')

clear('name','newStr1','path') 

%%
