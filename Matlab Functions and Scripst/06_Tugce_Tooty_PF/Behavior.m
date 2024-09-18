%% Behavior Analyse

% Variable -> data.bevaviour or data.bevaviour_bin (Averaging time epochs at desired sample rate)
% - first  cell row
%   .full record.

% - second cell row - from full session
%   . Row 1: index values for each event
%   . Row 2: time in samples for each event
%   . Row 3: percentage for each event

% - third cell row
%   . Column 1: Baseline
%   . Column 2: CS-Trial | Inter trial periods (ITI-epoch)

% - fourth cell line
%   . Imobillity lower threshold. Binary: 1 if <= threshold & 0 if > threshold

% - fifth cell line
%   . Row 1: index values for each event
%   . Row 2: time in samples for each event
%   . Row 3: percentage for each event

% - sixth cell line
%   . Total time in sec

% - seventh cell line
%   . Percentage %


% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  12/2023
% Last update: 12/2023

%% Data
data = [];
data.behavior{1,1} = raw';

clear('raw');

%% Trial epochs - baseline and CS sound period
% Original sample rate
parameters.original_srate = 30;

% Freezing settings
parameters.thr_1 = 200; % lower threshold in percentage.
parameters.thr_2 = 1; % higher threshold in sec to consider freezing.


% 3) Habituation

parameters.baseline_time   = 180 * parameters.original_srate; %samples
parameters.number_CS_tones = 5;
parameters.CS_tone_time    = 10 * parameters.original_srate; %samples
parameters.number_ITI      = 5;
parameters.ITI_time        = 30 * parameters.original_srate; %samples

parameters.total_time      = parameters.baseline_time  + (parameters.number_CS_tones * parameters.CS_tone_time) + (parameters.number_ITI * parameters.ITI_time);


% 3) Extinction

% parameters.baseline_time   = 180 * parameters.original_srate; %samples
% parameters.number_CS_tones = 45;
% parameters.CS_tone_time    = 10 * parameters.original_srate; %samples
% parameters.number_ITI      = 45;
% parameters.ITI_time        = 30 * parameters.original_srate; %samples
%
% parameters.total_time      = parameters.baseline_time  + (parameters.number_CS_tones * parameters.CS_tone_time) + (parameters.number_ITI * parameters.ITI_time);

%% data.events timestamps

data.events = [];

% Baseline session/Extinction
idx_start = (0:parameters.CS_tone_time + parameters.ITI_time:(parameters.CS_tone_time + parameters.ITI_time)*parameters.number_CS_tones) + parameters.baseline_time ;
idx_start(end) = [];

% CS_Tones
data.events{1,1}(:,1) = idx_start';                                              % Begin in samples
data.events{1,1}(:,2) = (data.events{1,1}(:,1) + parameters.CS_tone_time) - 1;   % End in samples

% ITI
data.events{1,2}(:,1) = (idx_start(2:end)' - parameters.ITI_time);               % Begin in samples
data.events{1,2}(:,2) = idx_start(2:end)' - 1;                                   % End in samples

data.events{1,2}(end+1,1) = data.events{1, 1}(end,2) + 1;                        % Begin the last ITI epoch
data.events{1,2}(end,2) = data.events{1, 1}(end,2) + parameters.ITI_time;        % End the last ITI epoch


clear('idx_start')


%% Full session Freezing

% idx for freezing epochs
[data.behavior{2,1}(1,:), data.behavior{2,1}(2,:), ~] = ZeroOnesCount(data.behavior{1,1} <= parameters.thr_1); % first row = start events / second = time in samples
idx_to_remove = data.behavior{2, 1}(2,:) < parameters.thr_2 * parameters.original_srate; % Exclude values based on timing threshold
data.behavior{2, 1}(:,idx_to_remove) = [];

if ~isempty(data.behavior{2, 1})

    % idx for non-freezing epochs
    % First lets find all freezing idx
    f_idx(:,1) = data.behavior{2, 1}(1,:);
    f_idx(:,2) = data.behavior{2, 1}(1,:) + data.behavior{2, 1}(2,:);

    %Then, lets find the intervals who represents non-freezing
    data.behavior{2,2}(:,1) = [1 f_idx(1,1)-1]; % from beggining until the first freezing event
    data.behavior{2,2}(1,2:length(f_idx)) = f_idx(1:end-1,2)+1;                      % first row = start events
    data.behavior{2,2}(2,2:length(f_idx)) = (f_idx(2:end,1)-1) - f_idx(1:end-1,2)-1; % second = time in samples
    data.behavior{2,2}(3,:) = data.behavior{2, 2}(2,:)./parameters.original_srate;   % thirt row = time in second

    idx_to_remove_1 = data.behavior{2, 2}(2,:) < parameters.thr_2 * parameters.original_srate; % Exclude values based on timing threshold
    data.behavior{2, 2}(:,idx_to_remove_1) = [];

else
    data.behavior{2, 1} = [];
end


clear('f_idx','idx_to_remove','idx_to_remove_1')

%% Trial epochs - Extinction, Retrieval - baseline, CS and ITI epochs

% baseline
data.behavior{3,1}  = data.behavior{1,1}(1,1:data.events{1, 1}(1,1)-1);

% CS sound and ITI period
for ii = 1:size(data.events{1, 1},1)

    epochs{1,ii}  = data.behavior{1,1}(1,data.events{1, 1}(ii,1):data.events{1, 1}(ii,2)); % CS-Trials
    epochs{2,ii}  = data.behavior{1,1}(1,data.events{1, 2}(ii,1):data.events{1, 2}(ii,2)); % ITI-Trials

end

% Reshaped CS and ITI trials in correct order and add to data behavior epochs
CS_ITI_Trials = reshape(epochs,1,[]);
data.behavior(3,2:length(CS_ITI_Trials)+1) = CS_ITI_Trials;


% CS-Trials and ITI Freezing
for ii = 1:size(data.behavior,2)

    data.behavior{4,ii} = data.behavior{3,ii} <= parameters.thr_1;

    [data.behavior{5,ii}(1,:), data.behavior{5,ii}(2,:), ~] = ZeroOnesCount(data.behavior{4,ii});

    idx_to_remove = data.behavior{5, ii}(2,:) < parameters.thr_2 * parameters.original_srate;
    data.behavior{5, ii}(:,idx_to_remove) = [];
    data.behavior{5, ii}(3,:) = data.behavior{5, ii}(2,:)./parameters.original_srate;

    idx_to_remove = [];

    data.behavior{6, ii} = sum(data.behavior{5, ii}(3,:));
    data.behavior{7, ii} = (data.behavior{6, ii}.*100) ./ (length(data.behavior{3, ii}) ./ parameters.original_srate);

end

clear('ii','epochs','CS_ITI_Trials','parameters.thr_1','parameters.thr_2','idx_to_remove')

%% Organizing freezing and non-freezing time-epochs considering the full sextion

if ~isempty(data.behavior{2,1})

    % Freezing idx
    for ii = 1:size(data.behavior{2,1},2)

        if data.behavior{2,1}(1,ii) < parameters.baseline_time; % Do not use baseline idx
            continue
        end

        data.events_behavior{1,1}(ii,1) = data.behavior{2,1}(1,ii);                            % Start freezing
        data.events_behavior{1,1}(ii,2) = data.behavior{2,1}(1,ii) + data.behavior{2,1}(2,ii); % End freezing

    end

    % delete zeros rows
    data.events_behavior{1, 1}( ~any(data.events_behavior{1, 1},2), : ) = [];  %rows

    % non-Freezing idx
    data.events_behavior{1,2}(:,1) = data.events_behavior{1,1}(1:end-1,2)+1;
    data.events_behavior{1,2}(:,2) = data.events_behavior{1,1}(2:end,1)-1;

    % Exclude time epochs < 1s
    idx_to_remove = (data.events_behavior{1,2}(:,2) - data.events_behavior{1,2}(:,1)) < 1*parameters.original_srate;
    data.events_behavior{1,2}(idx_to_remove,:) = [];


    clear('ii','idx_to_remove')

    % Set index of freezing and non-freezing for baseline, CS-Tones or ITI events

    % data.events_behavior organization
    %   - Row 1: Baseline freezing
    %   - Row 2: Baseline non-freezing
    %   - Row 3: CS-TONE freezing
    %   - Row 4: CS-TONE not freezing
    %   - Row 5: ITI freezing
    %   - Row 6: ITI non-freezing

    %   - Columns trials

    %       Within each cell:
    %           - Rows: events
    %           - Column 1: start
    %           - Column 2: end



    % Baseline
    % Freezing idx
    baseline_idx_f = data.behavior{2,1}(1,:) < data.events{1, 1}(1,1);  % logical idx
    data.events_behavior{2,1}(:,1) = data.behavior{2,1}(1,baseline_idx_f);                                            % Start idx
    data.events_behavior{2,1}(:,2) = data.behavior{2,1}(1,baseline_idx_f)'  + data.behavior{2,1}(2,baseline_idx_f)';  % End idx

    % Non-freezing
    baseline_idx_nf = data.behavior{2,2}(1,:) < data.events{1, 1}(1,1); % logical idx
    data.events_behavior{3,1}(:,1) = data.behavior{2,2}(1,baseline_idx_nf);                                               % Start idx
    data.events_behavior{3,1}(:,2) = data.behavior{2,2}(1,baseline_idx_nf)' + data.behavior{2,2}(2,baseline_idx_nf)' - 1; % End idx


    % CS-Tones
    for ii = 1:size(data.events{1, 1},1)

        % CS-Tone
        % Freezing idx
        CS_idx = data.behavior{2,1}(1,:) >= data.events{1, 1}(ii,1) & data.behavior{2,1}(1,:) <= data.events{1, 1}(ii,2);   % logical idx inside CS-TONE
        data.events_behavior{4,ii}(:,1) = data.behavior{2,1}(1,CS_idx);                                       % Start idx
        data.events_behavior{4,ii}(:,2) = data.behavior{2,1}(1,CS_idx) + data.behavior{2,1}(2,CS_idx) - 1;    % End idx
        % Non-freezing
        CS_idx_1 = data.behavior{2,2}(1,:) >= data.events{1, 1}(ii,1) & data.behavior{2,2}(1,:) <= data.events{1, 1}(ii,2); % logical idx inside CS-TONE
        data.events_behavior{5,ii}(:,1) = data.behavior{2,2}(1,CS_idx_1);                                      % Start idx
        data.events_behavior{5,ii}(:,2) = data.behavior{2,2}(1,CS_idx_1) + data.behavior{2,2}(2,CS_idx_1) - 1; % End idx

    end


    % ITI
    for ii = 1:size(data.events{1, 2},1)

        % ITI
        % Freezing idx
        ITI_idx = data.behavior{2,1}(1,:) >= data.events{1, 2}(ii,1) & data.behavior{2,1}(1,:) <= data.events{1, 2}(ii,2); % logical idx inside ITI
        data.events_behavior{6,ii}(:,1) = data.behavior{2,1}(1,ITI_idx);                                     % Start idx
        data.events_behavior{6,ii}(:,2) = data.behavior{2,1}(1,ITI_idx) + data.behavior{2,1}(2,ITI_idx) - 1; % End idx
        % Non-freezing
        ITI_idx = data.behavior{2,2}(1,:) >= data.events{1, 2}(ii,1) & data.behavior{2,2}(1,:) <= data.events{1, 2}(ii,2); % logical idx inside ITI
        data.events_behavior{7,ii}(:,1) = data.behavior{2,2}(1,ITI_idx);                                     % Start idx
        data.events_behavior{7,ii}(:,2) = data.behavior{2,2}(1,ITI_idx) + data.behavior{2,2}(2,ITI_idx) - 1; % End idx

    end

    clear('baseline_idx_f','baseline_idx_nf','CS_idx','CS_idx_1','ITI_idx','ii')

else

    data.events_behavior = [];

end

%% Temporary code used to save ITI timestamps
% It seems to be a redundant session, but it was the quickest way I found to group the variables (samples, seconds, and which trial happened). I need to adapt the code above.

if ~isempty(data.behavior{2,1})

    for ii = 1:size(data.events{1, 2},1)

        % ITI
        % Freezing idx
        ITI_idx = data.behavior{2,1}(1,:) >= data.events{1, 2}(ii,1) & data.behavior{2,1}(1,:) <= data.events{1, 2}(ii,2);            % logical idx inside ITI
        ITI_freezing_idx{1,ii}(:,1) = data.behavior{2,1}(1,ITI_idx);                                                                  % Start idx - samples
        ITI_freezing_idx{1,ii}(:,2) = data.behavior{2,1}(1,ITI_idx) + data.behavior{2,1}(2,ITI_idx) - 1;                              % End idx   - samples
        ITI_freezing_idx{1,ii}(:,3) = data.behavior{2,1}(1,ITI_idx)./parameters.original_srate;                                       % Start idx - seconds
        ITI_freezing_idx{1,ii}(:,4) = (data.behavior{2,1}(1,ITI_idx) + data.behavior{2,1}(2,ITI_idx) - 1)./parameters.original_srate; % Start idx - seconds
        ITI_freezing_idx{1,ii}(:,5) = ii;

    end


    temp1 = cat(1,ITI_freezing_idx{1,:});
    temp2 = repmat('Freezing',size(temp1,1),1); % For the FP spreadsheet , events labels need to be characters.

    % The FP recordings starts with a different time than zero. The time_factor variable simply corrects the times.
    % Each value was checked individually for each animal.
    % The timestamps here are being prepared for pMAT analysis.
    % Despite the standard organization between index Onset and index Offset, pMAT considers only the onset to organize the time windows to do PETH.
    % To analyze considering the offset, you need to change and put the offset
    % in the onset column and save it again. Okay, it seems very weird. But it's working like this.....( 3 -> onset / 4 -> offset)

    temp3 = temp1(:,[3 4]) + time_factor(ms,1); % save onset and offset
    temp4 = temp1(:,[4 4]) + time_factor(ms,1); % save offset and offset

    clear ('ii')

else

    temp2  = 0;
    temp3  = [ 0 0 ];
    temp4  = [ 0 0 ];

end

%% save spreedsheet with ITIs Timestamps to pMat

variables2save = {'Event','Onset','Offset'};
spreadsheet1 = table(string(temp2),temp3(:,1),temp3(:,2),'VariableNames',variables2save');
spreadsheet2 = table(string(temp2),temp4(:,1),temp4(:,2),'VariableNames',variables2save');

% Save *.xls

fprintf('\n Saving spreadsheets ...\n');

newStr = id(1:end-4);
path = '/Users/flavio/Desktop';
name1 = strcat(path,'/',newStr,'ITIevents_onset','.csv');
name2 = strcat(path,'/',newStr,'ITIevents_offset','.csv');

writetable(spreadsheet1,name1);
writetable(spreadsheet2,name2);

clear('name1','name2','newStr','path','temp1','temp2','temp3','temp4','spreadsheet1','spreadsheet2','variables2save','ITI_idx','ii','ITI_freezing_idx')

%% Temporary code used to save CS timestamps
% It seems to be a redundant session, but it was the quickest way I found to group the variables (samples, seconds, and trials). I need to adapt the code above.

if ~isempty(data.behavior{2,1})

    for ii = 1:size(data.events{1, 1},1)

        % CS
        % Freezing idx
        CS_idx = data.behavior{2,1}(1,:) >= data.events{1, 1}(ii,1) & data.behavior{2,1}(1,:) <= data.events{1, 1}(ii,2);           % logical idx inside CS
        CS_freezing_idx{1,ii}(:,1) = data.behavior{2,1}(1,CS_idx);                                                                  % Start idx - samples
        CS_freezing_idx{1,ii}(:,2) = data.behavior{2,1}(1,CS_idx) + data.behavior{2,1}(2,CS_idx) - 1;                               % End idx   - samples


        CS_freezing_idx{1,ii}(:,3) = data.behavior{2,1}(1,CS_idx)./parameters.original_srate;                                       % Start idx - seconds
        CS_freezing_idx{1,ii}(:,4) = (data.behavior{2,1}(1,CS_idx) + data.behavior{2,1}(2,CS_idx) - 1)./parameters.original_srate;  % Start idx - seconds
        CS_freezing_idx{1,ii}(:,5) = ii;

    end



    temp1 = cat(1,CS_freezing_idx{1,:});
    temp2 = repmat('Freezing',size(temp1,1),1); % For the FP spreadsheet , events labels need to be characters.

    % The FP recordings starts with a different time than zero. The time_factor variable simply corrects the times.
    % Each value was checked individually for each animal.
    % The timestamps here are being prepared for pMAT analysis.
    % Despite the standard organization between index Onset and index Offset, pMAT considers only the onset to organize the time windows to do PETH.
    % To analyze considering the offset, you need to change and put the offset
    % in the onset column and save it again. Okay, it seems very weird. But it's working like this.....( 3 -> onset / 4 -> offset)

    temp3 = temp1(:,[3 4]) + time_factor(ms,1); % save onset and offset
    temp4 = temp1(:,[4 4]) + time_factor(ms,1); % save offset and offset

else

    temp2  = 0;
    temp3  = [ 0 0 ];
    temp4  = [ 0 0 ];

end

%% save spreedsheet with CS Timestamps to pMat

variables2save = {'Event','Onset','Offset'};
spreadsheet1 = table(string(temp2),temp3(:,1),temp3(:,2),'VariableNames',variables2save');
spreadsheet2 = table(string(temp2),temp4(:,1),temp4(:,2),'VariableNames',variables2save');

% Save *.xls

fprintf('\n Saving spreadsheets ...\n');

newStr = id(1:end-4);
path = '/Users/flavio/Desktop';
name1 = strcat(path,'/',newStr,'CSevents_onset','.csv');
name2 = strcat(path,'/',newStr,'CSevents_offset','.csv');

writetable(spreadsheet1,name1);
writetable(spreadsheet2,name2);

clear('name1','name2','newStr','path','temp1','temp2','temp3','temp4','spreadsheet1','spreadsheet2','variables2save','CS_idx','ii','CS_freezing_idx')

%% Trial epochs - Binarized data in time epochs
% % Binarized data in time epochs just as analyzed in the laboratory.
%
% % -------------------------------------------
% % Important Note
%
% % At Plexon, Tugce delete some samples (before and after recording)...to fit? I don`t know why...
% % and exports data with bins of 200 ms (averaging 200 ms time epoch from their original sampling).
% % This represents an output similar to the conditioning box (5Hz sampling rate).
% % The data in its original sampling frequency recorded in Plexon versus "Binarized" by averaging,
% % may express slightly differently.
%
% % -------------------------------------------
%
% % Averaging time epochs at each 200ms(5hz)
% parameters.time_behav_bins = .2; %sec
%
% data_bins = 1:parameters.time_behav_bins * parameters.original_srate:length(data.behavior{1,1});
%
% for ii = 1:length(data_bins)-1
%     data.behavior_bins{1,1}(1,ii) = mean(data.behavior{1,1}(1,data_bins(ii):data_bins(ii+1)-1));
% end
%
%
%
% % Freezing settings
% parameters.thr_1 = 5; % lower threshold in percentage.
% parameters.thr_2 = 1; % higher threshold in sec to consider freezing.
%
% % Full session Freezing
%
% [data.behavior_bins{2,1}(1,:), data.behavior_bins{2,1}(2,:), ~] = ZeroOnesCount(data.behavior_bins{1,1} <= parameters.thr_1);
%
% idx_to_remove = data.behavior_bins{2, 1}(2,:) < (1/parameters.time_behav_bins);
% data.behavior_bins{2, 1}(:,idx_to_remove) = [];
% data.behavior_bins{2, 1}(3,:) = data.behavior_bins{2, 1}(2,:) ./ (1/parameters.time_behav_bins);
%
% %% Trial epochs - Extinction, retrieval - baseline, CS and ITI epochs
%
% % baseline
% data.behavior_bins{3,1}  = data.behavior_bins{1,1}(1,1:(data.events{2, 1}(1,1)-1)/(parameters.time_behav_bins * parameters.original_srate));
%
% % CS sound and ITI period
% for ii = 1:size(data.events{2, 1},1)
%
%         epochs_bins{1,ii}  = data.behavior_bins{1,1}(1,data.events{2, 1}(ii,1)/(parameters.time_behav_bins * parameters.original_srate) : (data.events{2, 1}(ii,2)-1)/(parameters.time_behav_bins * parameters.original_srate)); % CS-Trials
%         epochs_bins{2,ii}  = data.behavior_bins{1,1}(1,data.events{2, 2}(ii,1)/(parameters.time_behav_bins * parameters.original_srate) : (data.events{2, 2}(ii,2)-1)/(parameters.time_behav_bins * parameters.original_srate)); % ITI-Trials
%
% end
%
% % Reshaped CS and ITI trials in correct order and add to data behavior epochs
% CS_ITI_Trials_bins = reshape(epochs_bins,1,[]);
% data.behavior_bins(3,2:length(CS_ITI_Trials_bins)+1) = CS_ITI_Trials_bins;
%
%
% % CS-Trials and ITI Freezing
% for ii = 1:size(data.behavior_bins,2)
%
%     data.behavior_bins{4,ii} = data.behavior_bins{3,ii} <= parameters.thr_1;
%
%     [data.behavior_bins{5,ii}(1,:), data.behavior_bins{5,ii}(2,:), ~] = ZeroOnesCount(data.behavior_bins{4,ii});
%
%     idx_to_remove = data.behavior_bins{5, ii}(2,:) < parameters.thr_2 * (1/parameters.time_behav_bins);
%     data.behavior_bins{5, ii}(:,idx_to_remove) = [];
%     data.behavior_bins{5, ii}(3,:) = data.behavior_bins{5, ii}(2,:) ./ (1/parameters.time_behav_bins);
%
%     idx_to_remove = [];
%
%     data.behavior_bins{6, ii} = sum(data.behavior_bins{5, ii}(3,:));
%     data.behavior_bins{7, ii} = (data.behavior_bins{6, ii}.*100) ./ (length(data.behavior_bins{3,ii}) ./ (1/parameters.time_behav_bins));
%
% end
%
% clear('parameters.time_behav_bins','data_bins','ii','epochs_bins','CS_ITI_Trials_bins','parameters.thr_1','parameters.thr_2','idx_to_remove')

%% Graph based on Maren`s paper - nature communication
% Data freezing  Following the laboratory method, averages were calculated for every 5 time windows.

%Choose data:
data_2_use = data.behavior(:,1:end);

data.behavior_CS_epochs = zeros(1,(round(length(data_2_use)/5))/2 +1);
data.behavior_CS_epochs(1,1) = cell2mat(data_2_use(7, 1)); % baseline
data.behavior_CS_epochs(1,2:end) = mean(reshape(cell2mat(data_2_use(7, 2:2:end)),5,[]),1); % Averaging CS-Trials

clear ("data_2_use")

%% Select data to plot

% Movement
data_2_plot_1     = data.behavior{1,1};

% Time vector
behav_time_v = linspace(1,length(data.behavior{1,1})./parameters.original_srate,length(data.behavior{1,1}));

% CS indexes
cs_trial       = data.events{1, 1};

% Freezing indexes

if ~isempty(data.behavior{2,1})

freezing_ = cat(1,data.events_behavior{4,:});
freezing_start_CS = freezing_(:,1);
freezing_end_CS   = freezing_(:,2);

freezing_ = cat(1,data.events_behavior{6,:});
freezing_start_ITI = freezing_(:,1);
freezing_end_ITI   = freezing_(:,2);

else
    freezing_start_CS = nan;
    freezing_end_CS = nan;
    freezing_start_ITI = nan;
    freezing_end_ITI = nan;
end

% Freezing Percentage
data_2_plot_2  = cell2mat(data.behavior(7,:));

% Plot

figure
set(gcf,'color','w');
sc = [1,1,960,1200];
set(gcf, 'Position', sc);

%subplot(2,2,[1 3])
hold all
%sgtitle('Habituation')
%sgtitle('Exposure')
sgtitle(strcat(id(1:end-4),' - Extinction 1'))
%sgtitle('Retrieval')
%sgtitle('Renewal')


%plot(behav_bins_time_v,data_2_plot,'Color',[0.3, 0.3, 0.3, 0.3]) % raw data
plot(behav_time_v, data_2_plot_1,'linew', 1, 'Color',[0.8, 0.8, 0.8])
plot(behav_time_v,ones(1,length(behav_time_v)).*200,'k--')

if ~isempty(data.behavior{2,1})
    plot([behav_time_v(freezing_start_CS);behav_time_v(freezing_end_CS)], [ones(1,length(freezing_start_CS)).*4750;ones(1,length(freezing_end_CS)).*4750],'k-','linew', 5,'Color',[.6, 0, 0,.8])
    plot([behav_time_v(freezing_start_ITI);behav_time_v(freezing_end_ITI)], [ones(1,length(freezing_start_ITI)).*4750;ones(1,length(freezing_end_ITI)).*4750],'k-','linew', 5,'Color',	[0, 0.4470, 0.7410, .5])
end

plot([behav_time_v(cs_trial(:,1));behav_time_v(cs_trial(:,2))], [ones(1,length(cs_trial(:,1))).*5000;ones(1,length(cs_trial(:,2))).*5000],'k-','linew', 20,'Color',[1, .4, .4])

xlabel('Time (s)','FontSize',12), ylabel('Movement (%)','FontSize',12)
xlim([behav_time_v(1)-10 behav_time_v(end)])
ylim([0 5250])
%
% legend('Movement','Threshold','Freezing onset during CS-Tones','Freezing onset during ITIs','NumColumns',5,'Location','southoutside')
% legend('boxoff')

% subplot(2,2,2)
%
% % CS
% plot(data_2_plot_2(2:2:end),'-o','linew', 1, 'Color',[0.6, 0.6, 0.6, 0.6],'MarkerEdgeColor',[.6, 0, 0])
% ylabel('Freezing (%)','FontSize',12)
% xticklabels({'','CS 5','CS 10','CS 15','CS 20','CS 25','CS 30','CS 35','CS 40','CS 45'}) % CS-Trials
% %xticklabels({,'CS 1','CS 2','CS 3','CS 4','CS 5'}) % CS-Trials
% %xticklabels({'Baseline','CS 1','ITI','CS 2','ITI','CS 3','ITI','CS 4','ITI','CS 5','end'})       % CS & ITI-Trials
% %xticklabels({'Baseline','Blk 1','Blk 2','Blk 3','Blk 4','Blk 5','Blk 6','Blk 7','Blk 8','Blk 9'}) % Average block
% %xtickangle(90)
% xlim([0.2 length(data_2_plot_2(2:2:end))+1])
% ylim([-5 105])
% box off
% title('CS-Tones')
%
% subplot(2,2,4)
%
% %ITI
% plot(data_2_plot_2(3:2:end),'-o','linew', 1, 'Color',[0.6, 0.6, 0.6, 0.6],'MarkerEdgeColor',[0, 0, 0])
% ylabel('Freezing (%)','FontSize',12)
% xticklabels({'','ITI 5','ITI 10','ITI 15','ITI 20','ITI 25','ITI 30','ITI 35','ITI 40','ITI 45'}) % CS-Trials
% %xticklabels({,'CS 1','CS 2','CS 3','CS 4','CS 5'}) % CS-Trials
% %xticklabels({'Baseline','CS 1','ITI','CS 2','ITI','CS 3','ITI','CS 4','ITI','CS 5','end'})       % CS & ITI-Trials
% %xticklabels({'Baseline','Blk 1','Blk 2','Blk 3','Blk 4','Blk 5','Blk 6','Blk 7','Blk 8','Blk 9'}) % Average block
% %xtickangle(90)
% xlim([0.2 length(data_2_plot_2(3:2:end))+1])
% ylim([-5 105])
% box off
% title('ITIs')

% Clear

clear('behav_bins_time_v','cs_trial','data_2_plot_1','data_2_plot_2','freezing_end','freezing_start','ITI_freezing_idx','ITI_freezing','freezing_','behav_time_v','sc')


%% Save

%newStr = regexprep(files.id.name,'.mat','_');
%newStr = files.id(ms).name(1:end-8);
newStr = id(1:end-4);

path = '/Users/flavio/Desktop';
%path = files.FilesLoaded{1, 1}.folder;

name = strcat(path,'/',newStr,'_Extinction_1');

% save data
save(name,'data','parameters','-v7.3')

% save figure
saveas(gcf,name,'png')

close all
clear('name','newStr','path')

%% last update 18/01/2024 - 18:41
%  listening:
