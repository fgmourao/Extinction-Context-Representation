%% Organize Granger Freezing and Non-freezing epochs

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

%% Load files

% Load data pre-prepocessing and granger freezing and non-freezing

%% Pair channels combination
% Row 1: PL  --> IL
% Row 2: IL  --> PL
% Row 3: PL  --> HPC
% Row 4: HPC --> PL
% Row 5: IL  --> HPC
% Row 6: HPC --> IL

% All Possibles  combinations
combinations_ = zeros(size(mvgc.data{1,1},1)*2,2);
combinations_(1:2:6,:) = flip(nchoosek(1:size(mvgc.data,1),2),2);
combinations_(2:2:6,:) = nchoosek(1:size(mvgc.data,1),2);

%% Specs Baseline, freezing and non freezing
% This for-loop needs improvement due to redundancy in the baseline extraction data

mvgc.stats.Spect_2_12Hz = [];

% Baseline
for jj = 1:size(combinations_,1)
    mvgc.stats.Spect_2_12Hz{1,1}(jj,:) = mvgc.F_spect{1,1}(combinations_(jj,1), combinations_(jj,2),:);% baseline
end



% Define events that occurs within desired trials
  number_of_trials = 1;

  % First events
  number_of_events_freezing     = sum(data.events_behavior{2, 1}(:,1) <= data.events{2, 2}(number_of_trials,2));     % Considering first ITI idx for simetric CS presentation and ITI
  number_of_events_non_freezing = sum(data.events_behavior{2, 2}(:,1) <= data.events{2, 2}(number_of_trials,2));     % Considering first ITI idx for simetric non-CS presentation and ITI

  % Last events
%   number_of_events_freezing = sum(data.events_behavior{2, 1}(:,1) >= data.events{2, 2}(end-number_of_trials,2));     % Considering last ITI idx for simetric number os CS presentation and ITI
%   number_of_events_non_freezing = sum(data.events_behavior{2, 2}(:,1) <= data.events{2, 2}(number_of_trials,2));     % Considering last ITI idx for simetric non-CS presentation and ITI


% Freezing within CS
if isempty(data.events_behavior{4, 1})
    mvgc.stats.Spect_2_12Hz{2,1} = zeros(size(combinations_,1),3,mvgc.parameters.fres);
else
    %for ii = 1:length(data.events_behavior{4, 1})                        % Loop for all events
    for ii = 1:sum(data.events_behavior{4, 1}<=number_of_events_freezing) % Loop for first or last events
        for jj = 1:size(combinations_,1)
            mvgc.stats.Spect_2_12Hz{2,1}(jj,ii,:) = mvgc.F_spect{2,data.events_behavior{4, 1}(ii)}(combinations_(jj,1), combinations_(jj,2),:);      % All or First events
            %mvgc.stats.Spect_2_12Hz{2,1}(jj,ii,:) = mvgc.F_spect{2,data.events_behavior{4, 1}(end-ii)}(combinations_(jj,1), combinations_(jj,2),:); % Last events
            
        end
    end
end

% Freezing outside CS
if isempty(data.events_behavior{5, 1})
    mvgc.stats.Spect_2_12Hz{2,2} = zeros(size(combinations_,1),3,mvgc.parameters.fres);
else
    %for ii = 1:length(data.events_behavior{5, 1})                        % Loop for all events
    for ii = 1:sum(data.events_behavior{5, 1}<=number_of_events_freezing) % Loop for first or last events
        for jj = 1:size(combinations_,1)
            mvgc.stats.Spect_2_12Hz{2,2}(jj,ii,:) = mvgc.F_spect{2,data.events_behavior{5, 1}(ii)}(combinations_(jj,1), combinations_(jj,2),:);      % All or First events
            %mvgc.stats.Spect_2_12Hz{2,2}(jj,ii,:) = mvgc.F_spect{2,data.events_behavior{5, 1}(end-ii)}(combinations_(jj,1), combinations_(jj,2),:); % Last events

        end
    end
end

% Non-Freezing within CS
if isempty(data.events_behavior{4, 2})
    mvgc.stats.Spect_2_12Hz{3,1} = zeros(size(combinations_,1),3,mvgc.parameters.fres);
else
    %for ii = 1:length(data.events_behavior{4, 2})                            % Loop for all events
    for ii = 1:sum(data.events_behavior{4, 2}<=number_of_events_non_freezing) % Loop for first or last events
        for jj = 1:size(combinations_,1)
            mvgc.stats.Spect_2_12Hz{3,1}(jj,ii,:) = mvgc.F_spect{3,data.events_behavior{4, 2}(ii)}(combinations_(jj,1), combinations_(jj,2),:);     % All or First events
            %mvgc.stats.Spect_2_12Hz{3,1}(jj,ii,:) =mvgc.F_spect{3,data.events_behavior{4, 2}(end-ii)}(combinations_(jj,1), combinations_(jj,2),:); % Last events
        end
    end
end

% Non-Freezing outside CS
if isempty(data.events_behavior{5, 2})
    mvgc.stats.Spect_2_12Hz{3,2} = zeros(size(combinations_,1),3,mvgc.parameters.fres);
else
    %for ii = 1:length(data.events_behavior{5, 2})                            % Loop for all events
    for ii = 1:sum(data.events_behavior{5, 2}<=number_of_events_non_freezing) % Loop for first or last events
        for jj = 1:size(combinations_,1)
            mvgc.stats.Spect_2_12Hz{3,2}(jj,ii,:) = mvgc.F_spect{3,data.events_behavior{5, 2}(ii)}(combinations_(jj,1), combinations_(jj,2),:);      % All or First events
            %mvgc.stats.Spect_2_12Hz{3,2}(jj,ii,:) = mvgc.F_spect{3,data.events_behavior{5, 2}(end-ii)}(combinations_(jj,1), combinations_(jj,2),:); % Last events
        end
    end
end

clear('ii','jj')

%% Mean trials

mvgc.stats.Spect_2_12Hz_mean = [];

% Mean trials - CS and ITI
mvgc.stats.Spect_2_12Hz_mean{1,1} = mvgc.stats.Spect_2_12Hz{1,1};                    % Baseline - Redundant but helped in the code below
mvgc.stats.Spect_2_12Hz_mean{2,1} = squeeze(mean(mvgc.stats.Spect_2_12Hz{2,1},2));   % Freezing within CS
mvgc.stats.Spect_2_12Hz_mean{2,2} = squeeze(mean(mvgc.stats.Spect_2_12Hz{2,2},2));   % Freezing outside CS
mvgc.stats.Spect_2_12Hz_mean{3,1} = squeeze(mean(mvgc.stats.Spect_2_12Hz{3,1},2));   % Non-Freezing within CS
mvgc.stats.Spect_2_12Hz_mean{3,2} = squeeze(mean(mvgc.stats.Spect_2_12Hz{3,2},2));   % Non-Freezing outside CS


%% Data to Stats on Prism

mvgc.stats.Fint_3_6Hz = [];
mvgc.stats.Fint_6_9Hz = [];
% Matrices Organized into a single cell. Each row inside a cell represents a combination and each column a trial

% Baseline
for jj = 1:size(combinations_,1)

    mvgc.stats.Fint_3_6Hz{1,1}(jj,1) = mvgc.Fint_3_6Hz{1,1}(combinations_(jj,1), combinations_(jj,2)); % baseline
    mvgc.stats.Fint_6_9Hz{1,1}(jj,1) = mvgc.Fint_6_9Hz{1,1}(combinations_(jj,1), combinations_(jj,2)); % baseline

end

% Freezing within CS
%for ii = 1:length(data.events_behavior{4, 1})                        % Loop for all events
for ii = 1:sum(data.events_behavior{4, 1}<=number_of_events_freezing) % Loop for first or last events

    for jj = 1:size(combinations_,1)

        mvgc.stats.Fint_3_6Hz{2,1}(jj,ii) = mvgc.Fint_3_6Hz{2,data.events_behavior{4, 1}(ii)}(combinations_(jj,1), combinations_(jj,2));       % All or First events
        mvgc.stats.Fint_6_9Hz{2,1}(jj,ii) = mvgc.Fint_6_9Hz{2,data.events_behavior{4, 1}(ii)}(combinations_(jj,1), combinations_(jj,2));       % All or First events
        
        %mvgc.stats.Fint_3_6Hz{2,1}(jj,ii) = mvgc.Fint_3_6Hz{2,data.events_behavior{4, 1}(end-ii)}(combinations_(jj,1), combinations_(jj,2));  % Last events
        %mvgc.stats.Fint_6_9Hz{2,1}(jj,ii) = mvgc.Fint_6_9Hz{2,data.events_behavior{4, 1}(end-ii)}(combinations_(jj,1), combinations_(jj,2));  % Last events

    end
end

% Freezing outside CS
%for ii = 1:length(data.events_behavior{5, 1})                        % Loop for all events
for ii = 1:sum(data.events_behavior{5, 1}<=number_of_events_freezing) % Loop for first or last events
    for jj = 1:size(combinations_,1)

        mvgc.stats.Fint_3_6Hz{2,2}(jj,ii) = mvgc.Fint_3_6Hz{2,data.events_behavior{5, 1}(ii)}(combinations_(jj,1), combinations_(jj,2));       % All or First events
        mvgc.stats.Fint_6_9Hz{2,2}(jj,ii) = mvgc.Fint_6_9Hz{2,data.events_behavior{5, 1}(ii)}(combinations_(jj,1), combinations_(jj,2));       % All or First events

        %mvgc.stats.Fint_3_6Hz{2,2}(jj,ii) = mvgc.Fint_3_6Hz{2,data.events_behavior{5, 1}(end-ii)}(combinations_(jj,1), combinations_(jj,2));  % Last events
        %mvgc.stats.Fint_6_9Hz{2,2}(jj,ii) = mvgc.Fint_6_9Hz{2,data.events_behavior{5, 1}(end-ii)}(combinations_(jj,1), combinations_(jj,2));  % Last events

    end
end


% Non-Freezing within CS
%for ii = 1:length(data.events_behavior{4, 2})                            % Loop for all events
for ii = 1:sum(data.events_behavior{4, 2}<=number_of_events_non_freezing) % Loop for first or last events
    for jj = 1:size(combinations_,1)

        mvgc.stats.Fint_3_6Hz{3,1}(jj,ii) = mvgc.Fint_3_6Hz{3,data.events_behavior{4, 2}(ii)}(combinations_(jj,1), combinations_(jj,2));       % All or First events
        mvgc.stats.Fint_6_9Hz{3,1}(jj,ii) = mvgc.Fint_6_9Hz{3,data.events_behavior{4, 2}(ii)}(combinations_(jj,1), combinations_(jj,2));       % All or First events

        %mvgc.stats.Fint_3_6Hz{3,1}(jj,ii) = mvgc.Fint_3_6Hz{3,data.events_behavior{4, 2}(end-ii)}(combinations_(jj,1), combinations_(jj,2));  % Last events
        %mvgc.stats.Fint_6_9Hz{3,1}(jj,ii) = mvgc.Fint_6_9Hz{3,data.events_behavior{4, 2}(end-ii)}(combinations_(jj,1), combinations_(jj,2));  % Last events

    end
end

% Non-Freezing outside CS
%for ii = 1:length(data.events_behavior{5, 2})                            % Loop for all events
for ii = 1:sum(data.events_behavior{5, 2}<=number_of_events_non_freezing) % Loop for first or last events
    for jj = 1:size(combinations_,1)

        mvgc.stats.Fint_3_6Hz{3,2}(jj,ii) = mvgc.Fint_3_6Hz{3,data.events_behavior{5, 2}(ii)}(combinations_(jj,1), combinations_(jj,2));       % All or First events
        mvgc.stats.Fint_6_9Hz{3,2}(jj,ii) = mvgc.Fint_6_9Hz{3,data.events_behavior{5, 2}(ii)}(combinations_(jj,1), combinations_(jj,2));       % All or First events
        
        %mvgc.stats.Fint_3_6Hz{3,2}(jj,ii) = mvgc.Fint_3_6Hz{3,data.events_behavior{5, 2}(end-ii)}(combinations_(jj,1), combinations_(jj,2));  % Last events
        %mvgc.stats.Fint_6_9Hz{3,2}(jj,ii) = mvgc.Fint_6_9Hz{3,data.events_behavior{5, 2}(end-ii)}(combinations_(jj,1), combinations_(jj,2));  % Last events

    end
end


%% Matrices Organized into a single variable. in each cell each row represents a combination and the each column average trials
% I could make 'fancy' into a single variable, beautiful and sexy, but I'm tired and will confuse myself in the future with this because I'm dyslexic

mvgc.stats.Fint_3_6Hz_mean = [];
mvgc.stats.Fint_6_9Hz_mean = [];

% 3-6 Hz
% Mean trials
mvgc.stats.Fint_3_6Hz_mean{1,1}(:,1) = mvgc.stats.Fint_3_6Hz{1,1};            % baseline - Redundant but helped in the code below
mvgc.stats.Fint_3_6Hz_mean{2,1}(:,1) = mean(mvgc.stats.Fint_3_6Hz{2,1},2);    % Freezing within CS
mvgc.stats.Fint_3_6Hz_mean{3,1}(:,1) = mean(mvgc.stats.Fint_3_6Hz{3,1},2);    % Freezing outside CS
mvgc.stats.Fint_3_6Hz_mean{2,2}(:,1) = mean(mvgc.stats.Fint_3_6Hz{2,2},2);    % Non-Freezing within CS
mvgc.stats.Fint_3_6Hz_mean{3,2}(:,1) = mean(mvgc.stats.Fint_3_6Hz{3,2},2);    % Non-Freezing outside CS

% 6-9 Hz
% Mean trials
mvgc.stats.Fint_6_9Hz_mean{1,1}(:,1) = mvgc.stats.Fint_6_9Hz{1,1};            % baseline - Redundant but helped in the code below
mvgc.stats.Fint_6_9Hz_mean{2,1}(:,1) = mean(mvgc.stats.Fint_6_9Hz{2,1},2);    % Freezing within CS
mvgc.stats.Fint_6_9Hz_mean{3,1}(:,1) = mean(mvgc.stats.Fint_6_9Hz{3,1},2);    % Freezing outside CS
mvgc.stats.Fint_6_9Hz_mean{2,2}(:,1) = mean(mvgc.stats.Fint_6_9Hz{2,2},2);    % Non-Freezing within CS
mvgc.stats.Fint_6_9Hz_mean{3,2}(:,1) = mean(mvgc.stats.Fint_6_9Hz{3,2},2);    % Non-Freezing outside CS


%% peak frequency in range

mvgc.stats.Fint_3_6Hz_peak = [];
mvgc.stats.Fint_6_9Hz_peak = [];
mvgc.stats.Fint_3_6Hz_peak_mean = [];
mvgc.stats.Fint_6_9Hz_peak_mean = [];

% 3-6 Hz
% Peak on each trial frequency range
mvgc.stats.Fint_3_6Hz_peak{1,1} = max(mvgc.stats.Spect_2_12Hz{1,1}(:,mvgc.parameters.frex_idx_3_6Hz),[],2);    % baseline - Redundant but helped in the code below
mvgc.stats.Fint_3_6Hz_peak{2,1} = max(mvgc.stats.Spect_2_12Hz{2,1}(:,:,mvgc.parameters.frex_idx_3_6Hz),[],3);  % Freezing within CS
mvgc.stats.Fint_3_6Hz_peak{3,1} = max(mvgc.stats.Spect_2_12Hz{3,1}(:,:,mvgc.parameters.frex_idx_3_6Hz),[],3);  % Freezing outside CS
mvgc.stats.Fint_3_6Hz_peak{2,2} = max(mvgc.stats.Spect_2_12Hz{2,2}(:,:,mvgc.parameters.frex_idx_3_6Hz),[],3);  % Non-Freezing within CS
mvgc.stats.Fint_3_6Hz_peak{3,2} = max(mvgc.stats.Spect_2_12Hz{3,2}(:,:,mvgc.parameters.frex_idx_3_6Hz),[],3);  % Non-Freezing outside CS

% 6-9 Hz
% Peak value on each trial frequency range
mvgc.stats.Fint_6_9Hz_peak{1,1} = max(mvgc.stats.Spect_2_12Hz{1,1}(:,mvgc.parameters.frex_idx_6_9Hz(21:40)),[],2); % baseline - Redundant but helped in the code below
mvgc.stats.Fint_6_9Hz_peak{2,1} = max(mvgc.stats.Spect_2_12Hz{2,1}(:,:,mvgc.parameters.frex_idx_6_9Hz(21:40)),[],3); % Freezing within CS
mvgc.stats.Fint_6_9Hz_peak{3,1} = max(mvgc.stats.Spect_2_12Hz{3,1}(:,:,mvgc.parameters.frex_idx_6_9Hz(21:40)),[],3); % Freezing outside CS
mvgc.stats.Fint_6_9Hz_peak{2,2} = max(mvgc.stats.Spect_2_12Hz{2,2}(:,:,mvgc.parameters.frex_idx_6_9Hz(21:40)),[],3); % Non-Freezing within CS
mvgc.stats.Fint_6_9Hz_peak{3,2} = max(mvgc.stats.Spect_2_12Hz{3,2}(:,:,mvgc.parameters.frex_idx_6_9Hz(21:40)),[],3); % Non-Freezing outside CS


% 3-6 Hz
% Peak value Mean frequency range
mvgc.stats.Fint_3_6Hz_peak_mean{1,1} = max(mvgc.stats.Spect_2_12Hz_mean{1,1}(:,mvgc.parameters.frex_idx_3_6Hz),[],2);    % baseline - Redundant but helped in the code below
mvgc.stats.Fint_3_6Hz_peak_mean{2,1} = max(mvgc.stats.Spect_2_12Hz_mean{2,1}(:,mvgc.parameters.frex_idx_3_6Hz),[],2);  % Freezing within CS
mvgc.stats.Fint_3_6Hz_peak_mean{3,1} = max(mvgc.stats.Spect_2_12Hz_mean{3,1}(:,mvgc.parameters.frex_idx_3_6Hz),[],2);  % Freezing outside CS
mvgc.stats.Fint_3_6Hz_peak_mean{2,2} = max(mvgc.stats.Spect_2_12Hz_mean{2,2}(:,mvgc.parameters.frex_idx_3_6Hz),[],2);  % Non-Freezing within CS
mvgc.stats.Fint_3_6Hz_peak_mean{3,2} = max(mvgc.stats.Spect_2_12Hz_mean{3,2}(:,mvgc.parameters.frex_idx_3_6Hz),[],2);  % Non-Freezing outside CS

% 6-9 Hz
% Peak value Mean frequency range
mvgc.stats.Fint_6_9Hz_peak_mean{1,1} = max(mvgc.stats.Spect_2_12Hz_mean{1,1}(:,mvgc.parameters.frex_idx_6_9Hz(21:40)),[],2); % baseline - Redundant but helped in the code below
mvgc.stats.Fint_6_9Hz_peak_mean{2,1} = max(mvgc.stats.Spect_2_12Hz_mean{2,1}(:,mvgc.parameters.frex_idx_6_9Hz(21:40)),[],2); % Freezing within CS
mvgc.stats.Fint_6_9Hz_peak_mean{3,1} = max(mvgc.stats.Spect_2_12Hz_mean{3,1}(:,mvgc.parameters.frex_idx_6_9Hz(21:40)),[],2); % Freezing outside CS
mvgc.stats.Fint_6_9Hz_peak_mean{2,2} = max(mvgc.stats.Spect_2_12Hz_mean{2,2}(:,mvgc.parameters.frex_idx_6_9Hz(21:40)),[],2); % Non-Freezing within CS
mvgc.stats.Fint_6_9Hz_peak_mean{3,2} = max(mvgc.stats.Spect_2_12Hz_mean{3,2}(:,mvgc.parameters.frex_idx_6_9Hz(21:40)),[],2); % Non-Freezing outside CS


%% Full width at half maximum
% for ii = 1:size(mvgc.stats.Spect_2_12Hz_2_plot_mean,1)
%     for jj = 1:size(mvgc.stats.Spect_2_12Hz_2_plot_mean,2)
% 
%         [mvgc.stats.FWHM_3_6Hz_samples(ii,jj), mvgc.stats.FWHM_3_6Hz_Hz(ii,jj)] = FWHM(squeeze(mvgc.stats.Spect_2_12Hz_2_plot_mean(ii,jj,mvgc.parameters.frex_idx_3_6Hz))',mvgc.parameters.frex_3_6Hz);
%         [mvgc.stats.FWHM_6_9Hz_samples(ii,jj), mvgc.stats.FWHM_6_9Hz_Hz(ii,jj)] = FWHM(squeeze(mvgc.stats.Spect_2_12Hz_2_plot_mean(ii,jj,mvgc.parameters.frex_idx_6_9Hz))',mvgc.parameters.frex_6_9Hz);
%    
%     end
% end


%% Save data

% Settings
%ms = 1;
Path    = '/Users/flavio/Desktop';
newStr1 = id(1:end-8);
name_ = strcat(Path,'/',newStr1,'_MVGC_Granger_freezing_epochs_5_first_CS');

% Save data
save(name_,'mvgc','-v7.3')

clear('name_','newStr1','path') 

%%
