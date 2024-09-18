%% Organize Power Weltch Freezing and Non-freezing epochs

% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  02/2024
% Last update: 02/2024

% Organizing the matrices into a single variable. Each row represents the mean of the chosen trials in a respective pair of channels."
% Respective plots. Spectra and graphs

%% Load files


%% Specs Baseline, freezing and non freezing
% This for-loop needs improvement due to redundancy in the baseline extraction data

%% define frequncy range

theta_range = [2 12];
pw.params.closestfreq_1  = find(theta_range(1)<pw.full_trial.freq_baseline & pw.full_trial.freq_baseline<theta_range(2));
freq1_1 = pw.full_trial.freq_baseline(pw.params.closestfreq_1,1);

range = [3, 6];
closestfreq_2  = find(range(1)<freq1_1 & freq1_1<range(2));
freq2_2 = freq1_1 (closestfreq_2,1);

range = [6, 9];
closestfreq_3  = find(range(1)<freq1_1 & freq1_1<range(2));
freq3_3 = freq1_1 (closestfreq_3,1);

%%
pw.stats_total_power_first_trials = [];
pw.stats_total_power_first_trials_3_6Hz =[];
pw.stats_total_power_first_trials_6_9Hz =[];
pw.stats_total_power_first_trials_3_6Hz_peak =[];
pw.stats_total_power_first_trials_6_9Hz_peak =[];

% Baseline
pw.stats_total_power_first_trials{1,1} = pw.stats_total_power{1,1};  % baseline
pw.stats_total_power_first_trials_3_6Hz{1,1} = mean(pw.stats_total_power{1,1}(:,closestfreq_2),2); 
pw.stats_total_power_first_trials_6_9Hz{1,1} = mean(pw.stats_total_power{1,1}(:,closestfreq_3),2);

pw.stats_total_power_first_trials_3_6Hz_peak{1,1} = max(pw.stats_total_power{1,1}(:,closestfreq_2),[],2); 
pw.stats_total_power_first_trials_6_9Hz_peak{1,1} = max(pw.stats_total_power{1,1}(:,closestfreq_3),[],2);

% Define events that occurs within desired trials
  number_of_trials = 5;

  % First events
  number_of_events_freezing     = sum(data.events_behavior{2, 1}(:,1) <= data.events{2, 2}(number_of_trials,2));     % Considering first ITI idx for simetric CS presentation and ITI
  number_of_events_non_freezing = sum(data.events_behavior{2, 2}(:,1) <= data.events{2, 2}(number_of_trials,2));     % Considering first ITI idx for simetric non-CS presentation and ITI

  % Last events
%   number_of_events_freezing = sum(data.events_behavior{2, 1}(:,1) >= data.events{2, 2}(end-number_of_trials,2));     % Considering last ITI idx for simetric number os CS presentation and ITI
%   number_of_events_non_freezing = sum(data.events_behavior{2, 2}(:,1) <= data.events{2, 2}(number_of_trials,2));     % Considering last ITI idx for simetric non-CS presentation and ITI


% Freezing within CS
if isempty(data.events_behavior{4, 1})
    pw.stats_total_power_first_trials{2,1} = zeros(3,length(pw.params.closestfreq_1),3);
else
    %for ii = 1:length(data.events_behavior{4, 1})                        % Loop for all events
    for ii = 1:sum(data.events_behavior{4, 1}<=number_of_events_freezing) % Loop for first or last events
            
        pw.stats_total_power_first_trials{2,1}(:,:,ii) = pw.stats_total_power{1,2}(:,:,data.events_behavior{4, 1}(ii));       % All or First events
        
        pw.stats_total_power_first_trials_3_6Hz{2,1}(:,:,ii) = mean(pw.stats_total_power{1,2}(:,closestfreq_2,data.events_behavior{4, 1}(ii)),2);       % All or First events
        pw.stats_total_power_first_trials_3_6Hz_peak{2,1}(:,:,ii) = max(pw.stats_total_power{1,2}(:,closestfreq_2,data.events_behavior{4, 1}(ii)),[],2);     % All or First events
        
        pw.stats_total_power_first_trials_6_9Hz{2,1}(:,:,ii) = mean(pw.stats_total_power{1,2}(:,closestfreq_3,data.events_behavior{4, 1}(ii)),2);       % All or First events
        pw.stats_total_power_first_trials_6_9Hz_peak{2,1}(:,:,ii) = max(pw.stats_total_power{1,2}(:,closestfreq_3,data.events_behavior{4, 1}(ii)),[],2);     % All or First events
              
        %pw.stats_total_power_first_trials{2,1}(:,:,ii) = pw.stats_total_power{1,1}(:,:,data.events_behavior{4, 1}(end-1ii)); % Last events
            
    end
end

% Freezing outside CS
if isempty(data.events_behavior{5, 1})
    pw.stats_total_power_first_trials{2,2} = zeros(3,length(pw.params.closestfreq_1),3);
else
    %for ii = 1:length(data.events_behavior{5, 1})                        % Loop for all events
    for ii = 1:sum(data.events_behavior{5, 1}<=number_of_events_freezing) % Loop for first or last events
        
        pw.stats_total_power_first_trials{2,2}(:,:,ii) = pw.stats_total_power{1,2}(:,:,data.events_behavior{5, 1}(ii));     % All or First events

        pw.stats_total_power_first_trials_3_6Hz{2,2}(:,:,ii) = mean(pw.stats_total_power{1,2}(:,closestfreq_2,data.events_behavior{5, 1}(ii)),2);       % All or First events
        pw.stats_total_power_first_trials_3_6Hz_peak{2,2}(:,:,ii) = max(pw.stats_total_power{1,2}(:,closestfreq_2,data.events_behavior{5, 1}(ii)),[],2);     % All or First events
        
        pw.stats_total_power_first_trials_6_9Hz{2,2}(:,:,ii) = mean(pw.stats_total_power{1,2}(:,closestfreq_3,data.events_behavior{5, 1}(ii)),2);       % All or First events
        pw.stats_total_power_first_trials_6_9Hz_peak{2,2}(:,:,ii) = max(pw.stats_total_power{1,2}(:,closestfreq_3,data.events_behavior{5, 1}(ii)),[],2);     % All or First events

        %pw.stats_total_power_first_trials{2,2}(:,:,ii) = pw.stats_total_power{1,3}(:,:,data.events_behavior{5, 1}(end-ii)); % Last events
    
    end
end

% Non-Freezing within CS
if isempty(data.events_behavior{4, 2})
    pw.stats_total_power_first_trials{3,1} = zeros(3,length(pw.params.closestfreq_1),3);
else
    %for ii = 1:length(data.events_behavior{4, 2})                            % Loop for all events
    for ii = 1:sum(data.events_behavior{4, 2}<=number_of_events_non_freezing) % Loop for first or last events

            pw.stats_total_power_first_trials{3,1}(:,:,ii) = pw.stats_total_power{1,3}(:,:,data.events_behavior{4, 2}(ii));     % All or First events
            
            pw.stats_total_power_first_trials_3_6Hz{3,1}(:,:,ii) = mean(pw.stats_total_power{1,3}(:,closestfreq_2,data.events_behavior{4, 2}(ii)),2);       % All or First events
            pw.stats_total_power_first_trials_3_6Hz_peak{3,1}(:,:,ii) = max(pw.stats_total_power{1,3}(:,closestfreq_2,data.events_behavior{4, 2}(ii)),[],2);     % All or First events
        
            pw.stats_total_power_first_trials_6_9Hz{3,1}(:,:,ii) = mean(pw.stats_total_power{1,3}(:,closestfreq_3,data.events_behavior{4, 2}(ii)),2);       % All or First events
            pw.stats_total_power_first_trials_6_9Hz_peak{3,1}(:,:,ii) = max(pw.stats_total_power{1,3}(:,closestfreq_3,data.events_behavior{4, 2}(ii)),[],2);     % All or First events            
           
            %pw.stats_total_power_first_trials{3,1}(:,:,ii) = pw.stats_total_power{1,3}(:,:,data.events_behavior{4, 2}(ii)); % Last events

    end
end

% Non-Freezing outside CS
if isempty(data.events_behavior{5, 2})
    pw.stats_total_power_first_trials{3,2} = zeros(3,length(pw.params.closestfreq_1),3);
else
    %for ii = 1:length(data.events_behavior{5, 2})                            % Loop for all events
    for ii = 1:sum(data.events_behavior{5, 2}<=number_of_events_non_freezing) % Loop for first or last events

            pw.stats_total_power_first_trials{3,2}(:,:,ii) = pw.stats_total_power{1,3}(:,:,data.events_behavior{5, 2}(ii));       % All or First events

            pw.stats_total_power_first_trials_3_6Hz{3,2}(:,:,ii) = mean(pw.stats_total_power{1,3}(:,closestfreq_2,data.events_behavior{5, 2}(ii)),2);       % All or First events
            pw.stats_total_power_first_trials_3_6Hz_peak{3,2}(:,:,ii) = max(pw.stats_total_power{1,3}(:,closestfreq_2,data.events_behavior{5, 2}(ii)),[],2);     % All or First events
        
            pw.stats_total_power_first_trials_6_9Hz{3,2}(:,:,ii) = mean(pw.stats_total_power{1,3}(:,closestfreq_3,data.events_behavior{5, 2}(ii)),2);       % All or First events
            pw.stats_total_power_first_trials_6_9Hz_peak{3,2}(:,:,ii) = max(pw.stats_total_power{1,3}(:,closestfreq_3,data.events_behavior{5, 2}(ii)),[],2);     % All or First events    

            %pw.stats_total_power_first_trials{3,2}(:,:,ii) = pw.stats_total_power{1,3}(:,:,data.events_behavior{5, 2}(ii)); % Last events

    end
end

clear('ii','jj')

%% Mean trials

pw.stats_total_power_first_trials_mean = [];
pw.stats_total_power_first_trials_mean_3_6Hz = [];
pw.stats_total_power_first_trials_mean_6_9Hz = [];

% Mean trials - CS and ITI
pw.stats_total_power_first_trials_mean{1,1} = pw.stats_total_power_first_trials{1,1};                    % Baseline - Redundant but helped in the code below
pw.stats_total_power_first_trials_mean{2,1} = squeeze(mean(pw.stats_total_power_first_trials{2,1},3));   % Freezing within CS
pw.stats_total_power_first_trials_mean{2,2} = squeeze(mean(pw.stats_total_power_first_trials{2,2},3));   % Freezing outside CS
pw.stats_total_power_first_trials_mean{3,1} = squeeze(mean(pw.stats_total_power_first_trials{3,1},3));   % Non-Freezing within CS
pw.stats_total_power_first_trials_mean{3,2} = squeeze(mean(pw.stats_total_power_first_trials{3,2},3));   % Non-Freezing outside CS

pw.stats_total_power_first_trials_mean_3_6Hz{1,1} = pw.stats_total_power_first_trials_3_6Hz{1,1};                    % Baseline - Redundant but helped in the code below
pw.stats_total_power_first_trials_mean_3_6Hz{2,1} = squeeze(mean(pw.stats_total_power_first_trials_3_6Hz{2,1},3));   % Freezing within CS
pw.stats_total_power_first_trials_mean_3_6Hz{2,2} = squeeze(mean(pw.stats_total_power_first_trials_3_6Hz{2,2},3));   % Freezing outside CS
pw.stats_total_power_first_trials_mean_3_6Hz{3,1} = squeeze(mean(pw.stats_total_power_first_trials_3_6Hz{3,1},3));   % Non-Freezing within CS
pw.stats_total_power_first_trials_mean_3_6Hz{3,2} = squeeze(mean(pw.stats_total_power_first_trials_3_6Hz{3,2},3));   % Non-Freezing outside CS

pw.stats_total_power_first_trials_mean_6_9Hz{1,1} = pw.stats_total_power_first_trials_6_9Hz{1,1};                    % Baseline - Redundant but helped in the code below
pw.stats_total_power_first_trials_mean_6_9Hz{2,1} = squeeze(mean(pw.stats_total_power_first_trials_6_9Hz{2,1},3));   % Freezing within CS
pw.stats_total_power_first_trials_mean_6_9Hz{2,2} = squeeze(mean(pw.stats_total_power_first_trials_6_9Hz{2,2},3));   % Freezing outside CS
pw.stats_total_power_first_trials_mean_6_9Hz{3,1} = squeeze(mean(pw.stats_total_power_first_trials_6_9Hz{3,1},3));   % Non-Freezing within CS
pw.stats_total_power_first_trials_mean_6_9Hz{3,2} = squeeze(mean(pw.stats_total_power_first_trials_6_9Hz{3,2},3));   % Non-Freezing outside CS

%% Mean trials peaks

pw.stats_total_power_first_trials_mean_3_6Hz_peak = [];
pw.stats_total_power_first_trials_mean_6_9Hz_peak = [];

pw.stats_total_power_first_trials_mean_3_6Hz_peak{1,1} = pw.stats_total_power_first_trials_3_6Hz_peak{1,1};                    % Baseline - Redundant but helped in the code below
pw.stats_total_power_first_trials_mean_3_6Hz_peak{2,1} = squeeze(mean(pw.stats_total_power_first_trials_3_6Hz_peak{2,1},3));   % Freezing within CS
pw.stats_total_power_first_trials_mean_3_6Hz_peak{2,2} = squeeze(mean(pw.stats_total_power_first_trials_3_6Hz_peak{2,2},3));   % Freezing outside CS
pw.stats_total_power_first_trials_mean_3_6Hz_peak{3,1} = squeeze(mean(pw.stats_total_power_first_trials_3_6Hz_peak{3,1},3));   % Non-Freezing within CS
pw.stats_total_power_first_trials_mean_3_6Hz_peak{3,2} = squeeze(mean(pw.stats_total_power_first_trials_3_6Hz_peak{3,2},3));   % Non-Freezing outside CS

pw.stats_total_power_first_trials_mean_6_9Hz_peak{1,1} = pw.stats_total_power_first_trials_6_9Hz_peak{1,1};                    % Baseline - Redundant but helped in the code below
pw.stats_total_power_first_trials_mean_6_9Hz_peak{2,1} = squeeze(mean(pw.stats_total_power_first_trials_6_9Hz_peak{2,1},3));   % Freezing within CS
pw.stats_total_power_first_trials_mean_6_9Hz_peak{2,2} = squeeze(mean(pw.stats_total_power_first_trials_6_9Hz_peak{2,2},3));   % Freezing outside CS
pw.stats_total_power_first_trials_mean_6_9Hz_peak{3,1} = squeeze(mean(pw.stats_total_power_first_trials_6_9Hz_peak{3,1},3));   % Non-Freezing within CS
pw.stats_total_power_first_trials_mean_6_9Hz_peak{3,2} = squeeze(mean(pw.stats_total_power_first_trials_6_9Hz_peak{3,2},3));   % Non-Freezing outside CS


%% Save data

% Settings
%ms = 1;
Path    = '/Users/flavio/Desktop';
newStr1 = id(1:end-8);
name_ = strcat(Path,'/',newStr1,'_pw_nFFT_32768_freezingg_epochs_5_first_CS');


% Save data
save(name_,'pw','-v7.3')

clear('name_','newStr1','path') 

%%
