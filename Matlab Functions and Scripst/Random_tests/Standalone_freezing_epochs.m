%% Standalone script organizes the freezing and non-freezing

% This standalone script organizes the freezing and non-freezing indices according to time criteria and whether they occurred within or outside the CS or ITI intervals.

% ----------------------------------
% IMPORTANT: All these ideas were developed after the initial preprocessing. Initially, the recording had been cut only for CS and ITI windows.
% It was then followed by processing the freezing and non-freezing windows as a whole.
% Now, a small piece of code (at line 120) has been attached to identify whether the freezing and non-freezing windows occurred within or outside of CS/ITI.
% ----------------------------------

% To reiterate:
% The processing of the recording was carried out for freezing and non-freezing windows without discriminating whether they were within or outside of CS/ITI.
% Until now, this curation must be done manually during the plotting and statistical analysis of the data.

% I chose to keep a standalone script until all time variables and decisions are made. In the future, these lines will be appended to the original scripts


% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  02/2024
% Last update:

%% Identifying the indices of freezing and non-freezing.

Baseline_time = data.events_behavior{1, 1}(1,1); % Defined by the fear protocol. all before first CS-Tone

data.events_behavior =[];

% Freezing idx
for ii = 1:size(data.behavior{2,1},2)

    if data.behavior{2,1}(1,ii) < Baseline_time * parameters.original_srate; % Do not use baseline idx
        continue
    end

    data.events_behavior{1,1}(ii,1) = data.behavior{2,1}(1,ii);                            % Start freezing
    data.events_behavior{1,1}(ii,2) = data.behavior{2,1}(1,ii) + data.behavior{2,1}(2,ii); % End freezezing

end

% delete zeros rows
data.events_behavior{1, 1}( ~any(data.events_behavior{1, 1},2), : ) = [];  %rows


% non-Freezing idx
data.events_behavior{1,2}(:,1) = data.events_behavior{1,1}(1:end-1,2)+1; % Start non-freezing
data.events_behavior{1,2}(:,2) = data.events_behavior{1,1}(2:end,1)-1;   % End non-freezing

% Exclude time epochs < 1s
idx_to_remove = (data.events_behavior{1,2}(:,2) - data.events_behavior{1,2}(:,1)) < 1*parameters.original_srate;
data.events_behavior{1,2}(idx_to_remove,:) = [];


%% Organizing the indices of freezing and non-freezing guided by the timing hypothesis.

% selected minimum time windows for each event
timewin = 3; % sec.

data.events_behavior{2, 1} = (data.events_behavior{1, 1}(:,2) - data.events_behavior{1, 1}(:,1))>=timewin.*parameters.decimated_srate;
data.events_behavior{2, 1} = data.events_behavior{1, 1}(data.events_behavior{2, 1},:);
data.events_behavior{2, 2} = (data.events_behavior{1, 2}(:,2) - data.events_behavior{1, 2}(:,1))>=timewin.*parameters.decimated_srate;
data.events_behavior{2, 2} = data.events_behavior{1, 2}(data.events_behavior{2, 2},:);

clear('ii','jj','timewin','data.events_behavior{2, 1}','data.events_behavior{2, 2}','time_win_pos','time_win_pos','time_win_pre' )

%% Organizing data - Considering only the freezing and non-freezing behavior


% Freezing onset
time_win_pre = 1.5; % sec
% Freezing offset
time_win_pos = 1.5; % sec


if ~isempty(data.lfp{9,1})
    data.lfp(9:12,:)  = [];
end


for ii = 1:size(data.lfp,2)
    for jj = 1:size(data.events_behavior{2, 1},1)

        % Freezing onset
        data.lfp{9,ii}(:,:,jj) = data.lfp{5,ii}(:,((data.events_behavior{2, 1}(jj,1) - time_win_pre * parameters.decimated_srate)+1 : (data.events_behavior{2, 1}(jj,1) + timewin * parameters.decimated_srate)));

        % Freezing offset
        if jj ~= size(data.events_behavior{2, 1},1) % ignoring last idx.
            data.lfp{10,ii}(:,:,jj) = data.lfp{5,ii}(:,((data.events_behavior{2, 1}(jj,2) - timewin * parameters.decimated_srate) : (data.events_behavior{2, 1}(jj,2) + time_win_pos * parameters.decimated_srate)-1));
        end

        % Freezing epochs
        data.lfp{11,ii}(:,:,jj) = data.lfp{5,ii}(:,(data.events_behavior{2, 1}(jj,1):data.events_behavior{2, 1}(jj,1) + timewin * parameters.decimated_srate - 1));


        if isempty(data.lfp{1,ii}) % The raw data has never been filtered until now.
            continue
        else
            data.lfp{3,ii}(:,:,jj) = data.lfp{1,ii}(:,data.events{2, 1}(jj,1):data.events{2, 1}(jj,2)-1);
        end
    end

    for jj = 1:size(data.events_behavior{2, 2},1)

        % non-Freezing epochs
        data.lfp{12,ii}(:,:,jj) = data.lfp{5,ii}(:,(data.events_behavior{2, 2}(jj,1):data.events_behavior{2, 2}(jj,1) + timewin * parameters.decimated_srate - 1));

        if isempty(data.lfp{1,ii}) % The raw data has never been filtered until now.
            continue
        else
            data.lfp{3,ii}(:,:,jj) = data.lfp{1,ii}(:,data.events{2, 1}(jj,1):data.events{2, 1}(jj,2)-1);
        end
    end

end

%% Identifying whether each freezing and non-freezing event occurred within the CS or within the ITI

% Searching for logical values. Freezing event started within the CS epoch or within the ITI epoch
for ii = 1:size(data.events_behavior{2, 1},1)

    data.events_behavior{3, 1}(ii,1) = sum(data.events_behavior{2, 1}(ii,1)>=data.events{2, 1}(:,1) & data.events_behavior{2, 1}(ii,1)<=data.events{2, 1}(:,2)); % logical values. Freezing event started within the CS epoch
    data.events_behavior{3, 1}(ii,2) = sum(data.events_behavior{2, 1}(ii,1)>=data.events{2, 2}(:,1) & data.events_behavior{2, 1}(ii,1)<=data.events{2, 2}(:,2)); % logical values. Freezing event started within the ITI epoch

end

% Event index of Freezing who start within the CS epoch or within the ITI epoch
if ~isempty(data.events_behavior{2, 1})
    data.events_behavior{4, 1}(:,1) = find(data.events_behavior{3, 1}(:,1)); % event index who start within the CS epoch
    data.events_behavior{5, 1}(:,1) = find(data.events_behavior{3, 1}(:,2)); % event index who start within the ITI epoch

end


% Searching for logical values. Non-Freezing event started within the CS epoch or within the ITI epoch
for ii = 1:size(data.events_behavior{2, 2},1)

    data.events_behavior{3, 2}(ii,1) = sum(data.events_behavior{2, 2}(ii,1)>=data.events{2, 1}(:,1) & data.events_behavior{2, 2}(ii,1)<=data.events{2, 1}(:,2)); % logical values. Non-freezing event started within the CS epoch
    data.events_behavior{3, 2}(ii,2) = sum(data.events_behavior{2, 2}(ii,1)>=data.events{2, 2}(:,1) & data.events_behavior{2, 2}(ii,1)<=data.events{2, 2}(:,2)); % logical values. Non-freezing event started within the ITI epoch

end

% Event index of Non-Freezing who start within the CS epoch or within the ITI epoch
if ~isempty(data.events_behavior{2, 2})
    data.events_behavior{4, 2}(:,1) = find(data.events_behavior{3, 2}(:,1)); % event index who start within the CS epoch
    data.events_behavior{5, 2}(:,1) = find(data.events_behavior{3, 2}(:,2)); % event index who start within the ITI epoch

end

%% last update 17/02/2024 -
%  listening: