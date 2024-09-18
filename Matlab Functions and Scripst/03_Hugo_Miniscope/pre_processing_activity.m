
%% Analysis of cellular activity. Calcium imaging experiment. Miniscope.
%  Pre-Processing - Sorting Session, events and cells

%  After preprocessing, the data.activity from both sessions_activity were merged and saved in the same spreadsheet.
%  the table includes time (seconds) where the event occurs, the cell ID that produced the event, and the amplitude (value) of each event.

%  The MedPC trigger starts the recording, so 0s on the miniscope recording is also 0s on the freezing output


% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  02/2024
% Last update: 03/2024

%% Experimental sessions_activity

% 1) Baseline -> 300

% 2) Consolidation -> 20 min

% 3) Extinction
%    - Baseline - 180s
%    - 45 CS-Tones -> 10 sec each
%    - 45 ITI -> 60 sec each

% 4) Retrieval
%    - Baseline - 180s
%    - 45 CS-Tones -> 10 sec each
%    - 45 ITI -> 60 sec each

%% Load and sort sessions_activity. 

% First: load spreedsheet
data_excel = Rat1cellstimespemps; % Spreedsheet name

% Let`s sorting by time to find the gap between sessions_activity and thus separate them
time_sort = sortrows(data_excel,1,'ascend');

% index where each session started. 
% Considering that each session lasts a minimum of 300 seconds, 250 was a reasonable number.
session_idx_ = find([1 ;diff(time_sort.Times) > 250]);

% column 1 --> extinction
% Column 2 --> retrieval

data = {};

for jj = 1:length(session_idx_)
    for ii = 1:size(time_sort,2)

        % Session

        if jj == length(session_idx_)
            data.activity{jj,1}(:,ii) = table2array(time_sort(session_idx_(jj):end,ii));

        else
            data.activity{jj,1}(:,ii) = table2array(time_sort(session_idx_(jj):session_idx_(jj+1)-1,ii));

        end
    end
end


% Correcting the time for each session
for ii = 1:size(data.activity,1)
    data.activity{ii,1}(:,1) = data.activity{ii,1}(:,1) - data.activity{ii,1}(1,1);
end


clear('session_idx_','ii','jj')

%% Settings from behavior experimental sessions_activity

parameters.sampling_rate = 5;

% 1) Baseline -> 300
parameters.baseline.total_time = 300; % samples
parameters.baseline.blocks     = 60;  % samples

% 2) Consolidation -> 20 min
parameters.consolidation.total_time = 1200; % samples
parameters.consolidation.blocks = 300;  % samples

% 3) Extinction

parameters.extinction.baseline_time   = 180; %samples
parameters.extinction.number_CS_tones = 45;
parameters.extinction.CS_tone_time    = 10; %samples
parameters.extinction.number_ITI      = 45;
parameters.extinction.ITI_time        = 30; %samples
parameters.extinction.total_time      = parameters.extinction.baseline_time  + (parameters.extinction.number_CS_tones * parameters.extinction.CS_tone_time) + (parameters.extinction.number_ITI * parameters.extinction.ITI_time);

% 3) retrieval
parameters.retrieval.baseline_time   = 180; %samples
parameters.retrieval.number_CS_tones = 45;
parameters.retrieval.CS_tone_time    = 10; %samples
parameters.retrieval.number_ITI      = 45;
parameters.retrieval.ITI_time        = 30; %samples
parameters.retrieval.total_time      = parameters.retrieval.baseline_time  + (parameters.retrieval.number_CS_tones * parameters.retrieval.CS_tone_time) + (parameters.retrieval.number_ITI * parameters.retrieval.ITI_time);


%% data.events timestamps

data.events = [];


% Baseline
data.events{1,1}(:,1) = 0:parameters.baseline.blocks:parameters.baseline.total_time;
data.events{1,1}(:,2) = data.events{1,1}(:,1) + (parameters.baseline.blocks)-1;
data.events{1,1}(end,:) = [];


% Consolidation
data.events{2,1}(:,1) = 0:parameters.consolidation.blocks:parameters.consolidation.total_time;
data.events{2,1}(:,2) = data.events{2,1}(:,1) + (parameters.consolidation.blocks)-1;
data.events{2,1}(end,:) = [];


% Extinction
idx_start = (0:parameters.extinction.CS_tone_time + parameters.extinction.ITI_time:(parameters.extinction.CS_tone_time + parameters.extinction.ITI_time)*45) + parameters.extinction.baseline_time ;
idx_start(end) = [];

% CS_Tones
data.events{3,1}(:,1) = idx_start';                                                         % Begin in samples
data.events{3,1}(:,2) = (data.events{3,1}(:,1) + parameters.extinction.CS_tone_time) - 1;   % End in samples

% ITI
data.events{3,2}(:,1) = (idx_start(2:end)' - parameters.extinction.ITI_time);               % Begin in samples
data.events{3,2}(:,2) = idx_start(2:end)' - 1;                                              % End in samples

data.events{3,2}(end+1,1) = data.events{3, 1}(end,2) + 1;                                   % Begin the last ITI epoch
data.events{3,2}(end,2) = data.events{3, 1}(end,2) + parameters.extinction.ITI_time;        % End the last ITI epoch


% Retrieval
% CS_Tones
data.events{4,1}(:,1) = idx_start';                                                         % Begin in samples
data.events{4,1}(:,2) = (data.events{4,1}(:,1) + parameters.extinction.CS_tone_time) - 1;   % End in samples

% ITI
data.events{4,2}(:,1) = (idx_start(2:end)' - parameters.extinction.ITI_time);               % Begin in samples
data.events{4,2}(:,2) = idx_start(2:end)' - 1;                                              % End in samples

data.events{4,2}(end+1,1) = data.events{4, 1}(end,2) + 1;                                   % Begin the last ITI epoch
data.events{4,2}(end,2) = data.events{4, 1}(end,2) + parameters.extinction.ITI_time;        % End the last ITI epoch

clear('idx_start')

%% Organizing data - Considering the baseline CS and ITI Events

% data.activity

% Row 1   : Baseline
% Column 1: Entire trial
% Column 2: Blocks with time determined by events

% Row 1   : Consolidation
% Column 1: Entire trial
% Column 2: Blocks with time determined by events

% Row 3   : Extinction
% Column 1: Entire trial
% Column 2: Baseline
% From Column 3 and Odd columns: CS-tone
% From Column 4 and Even columns: ITI

% Row 4: Retrieval
% Column 1: Entire trial
% Column 2: Baseline
% From Column 3 and Odd columns: CS-tone
% From Column 4 and Even columns: ITI

%data.activity = [];

% Baseline
for ii = 1:size(data.events{1, 1},1)
    data.activity{1,ii+1} = data.activity{1,1}(data.activity{1,1}(:,1) >= data.events{1, 1}(ii,1) & data.activity{1,1}(:,1) <= data.events{1, 1}(ii,2),:);
    
    if ii > 1
        data.activity{1,ii+1}(:,1) = data.activity{1,ii+1}(:,1) - data.events{1, 1}(ii,1);
    end
    
end


% Consolidation
for ii = 1:size(data.events{2, 1},1)
    data.activity{2,ii+1} = data.activity{2,1}(data.activity{2,1}(:,1) >= data.events{2, 1}(ii,1) & data.activity{2,1}(:,1) <= data.events{2, 1}(ii,2),:);
    
    if ii > 1
        data.activity{2,ii+1}(:,1) = data.activity{2,ii+1}(:,1) - data.events{2, 1}(ii,1);
    end

end


% Extinction
% Baseline
data.activity{3,2} = data.activity{3,1}(data.activity{3,1}(:,1) < parameters.extinction.baseline_time,:);
% CS-TONES
for ii = (1:2:(size(data.events{3, 1},1) + size(data.events{3, 2},1)))+2
    data.activity{3,ii} = data.activity{3,1} (data.activity{3,1}(:,1) >= data.events{3, 1}(floor(ii/2),1) & data.activity{3,1}(:,1) <= data.events{3, 1}(floor(ii/2),2),:);
    data.activity{3,ii+1} = data.activity{3,1}(data.activity{3,1}(:,1) >= data.events{3, 2}(floor(ii/2),1) & data.activity{3,1}(:,1) <= data.events{3, 2}(floor(ii/2),2),:);
end

% Correcting time for each trial
% CS-Tone
for ii = 3:2:(parameters.extinction.number_CS_tones+parameters.extinction.number_CS_tones)+2
        data.activity{3,ii}(:,1) = data.activity{3,ii}(:,1) - data.events{3, 1}(floor(ii/2),1);
end

% ITI
for ii = 4:2:(parameters.extinction.number_ITI+parameters.extinction.number_ITI)+2
        data.activity{3,ii}(:,1) = data.activity{3,ii}(:,1) - data.events{3, 2}((ii/2)-1,1);
end


% Retrieval
% Baseline
data.activity{4,2} = data.activity{4,1}(data.activity{4,1}(:,1) < parameters.retrieval.baseline_time,:);
% CS-TONES
for ii = (1:2:(size(data.events{4, 1},1) + size(data.events{4, 2},1)))+2
    data.activity{4,ii} = data.activity{4,1} (data.activity{4,1}(:,1) >= data.events{4, 1}(floor(ii/2),1) & data.activity{4,1}(:,1) <= data.events{4, 1}(floor(ii/2),2),:);
    data.activity{4,ii+1} = data.activity{4,1}(data.activity{4,1}(:,1) >= data.events{4, 2}(floor(ii/2),1) & data.activity{4,1}(:,1) <= data.events{4, 2}(floor(ii/2),2),:);

end

% Correcting time for each trial
% CS-Tone
for ii = 3:2:(parameters.extinction.number_CS_tones+parameters.extinction.number_CS_tones)+2
        data.activity{4,ii}(:,1) = data.activity{4,ii}(:,1) - data.events{4, 1}(floor(ii/2),1);
end

% ITI
for ii = 4:2:(parameters.extinction.number_ITI+parameters.extinction.number_ITI)+2
        data.activity{4,ii}(:,1) = data.activity{4,ii}(:,1) - data.events{4, 2}((ii/2)-1,1);
end

clear('ii')

%% Baseline - Sorting Cells activity and timestamps

% Categorizing the ID of each cell and the number of firings (events) of each one.
% GC --> numer of events
% GR --> cell ID

data.activity_sorted_baseline = [];

for ii = 1:size(data.events{1, 1},1)+1
    [data.activity_sorted_baseline{1,ii}(:,2),data.activity_sorted_baseline{1,ii}(:,1)] = groupcounts(data.activity{1,ii}(:,2));
end

% In Each cell row
% -->  row: cells
% -->  columns: time

for ii = 1:size(data.events{1, 1},1)+1
    data.activity_sorted_baseline{2,ii} = NaN(max(data.activity_sorted_baseline{1,1}(:,1)),max(data.activity_sorted_baseline{1,ii}(:,2))); % --> time in samples
    data.activity_sorted_baseline{3,ii} = NaN(max(data.activity_sorted_baseline{1,1}(:,1)),max(data.activity_sorted_baseline{1,ii}(:,2))); % --> time in seconds
    data.activity_sorted_baseline{4,ii} = NaN(max(data.activity_sorted_baseline{1,1}(:,1)),max(data.activity_sorted_baseline{1,ii}(:,2))); % --> amplitude events

    if ii == 1
        data.activity_sorted_baseline{5,ii} = zeros(size(data.activity_sorted_baseline{1, 1},1),parameters.baseline.total_time * parameters.sampling_rate); % --> entire session binary format
        data.activity_sorted_baseline{6,ii} = zeros(size(data.activity_sorted_baseline{1, 1},1),parameters.baseline.total_time * parameters.sampling_rate); % --> entire session amplitude values
    elseif ii>1
        data.activity_sorted_baseline{5,ii} = zeros(size(data.activity_sorted_baseline{1, 1},1),parameters.baseline.blocks * parameters.sampling_rate); % --> entire session binary format
        data.activity_sorted_baseline{6,ii} = zeros(size(data.activity_sorted_baseline{1, 1},1),parameters.baseline.blocks * parameters.sampling_rate); % --> entire session amplitude values
    end

end

for jj = 1:size(data.events{1, 1},1)+1
    for ii = 1:size(data.activity_sorted_baseline{1,jj},1)

        temp_data = data.activity{1,jj}(data.activity{1,jj}(:,2) == ii, 1)';

        % time in samples
        data.activity_sorted_baseline{2,jj}(data.activity_sorted_baseline{1, jj}(ii,1),1:length(temp_data)) = round(temp_data * parameters.sampling_rate);
        data.activity_sorted_baseline{2,jj}(data.activity_sorted_baseline{1, jj}(ii,1),data.activity{2,1}(ii,:) == 0) = 1; % As the data.activity was transformed from seconds to samples, times of 0 were converted to sample 1.

        % Time in seconds
        data.activity_sorted_baseline{3,jj}(data.activity_sorted_baseline{1, jj}(ii,1),1:length(temp_data)) = temp_data;

        % Binary events over time
        temp_data((temp_data * parameters.sampling_rate) == 0) = 1;     % As the data.activity was transformed from seconds to samples, times of 0 were converted to sample 1.
        data.activity_sorted_baseline{5,jj}(data.activity_sorted_baseline{1, jj}(ii,1),ceil(temp_data * parameters.sampling_rate)) = 1;  % time sample to binary

        %  Amplitude of the event
        temp_data_amplitude = data.activity{1,jj}(data.activity{1,jj}(:,2) == ii, 3)';
        data.activity_sorted_baseline{4,jj}(data.activity_sorted_baseline{1, jj}(ii,1),1:length(temp_data_amplitude)) = temp_data_amplitude;

        %  Amplitude of the event over time
        data.activity_sorted_baseline{6,jj}(data.activity_sorted_baseline{1, jj}(ii,1),ceil(temp_data * parameters.sampling_rate)) = temp_data_amplitude;  % time sample to binary

    end

    % Time vector in seconds full session
    data.activity_sorted_baseline{7,jj} = linspace(0,length(data.activity_sorted_consolidation{6,jj})/parameters.sampling_rate,length(data.activity_sorted_baseline{6,jj}));

end

clear('temp_data','temp_data_amplitude','ii','jj')

%% Consolidation - Sorting Cells activity and timestamps

% Categorizing the ID of each cell and the number of firings (events) of each one.
% GC --> numer of events
% GR --> cell ID

data.activity_sorted_consolidation = [];

for ii = 1:size(data.events{2, 1},1)+1
    [data.activity_sorted_consolidation{1,ii}(:,2),data.activity_sorted_consolidation{1,ii}(:,1)] = groupcounts(data.activity{2,ii}(:,2));
end

% In Each cell row
% -->  row: cells
% -->  columns: time

for ii = 1:size(data.events{2, 1},1)+1
    data.activity_sorted_consolidation{2,ii} = NaN(max(data.activity_sorted_consolidation{1,1}(:,1)),max(data.activity_sorted_consolidation{1,ii}(:,2))); % --> time in samples
    data.activity_sorted_consolidation{3,ii} = NaN(max(data.activity_sorted_consolidation{1,1}(:,1)),max(data.activity_sorted_consolidation{1,ii}(:,2))); % --> time in seconds
    data.activity_sorted_consolidation{4,ii} = NaN(max(data.activity_sorted_consolidation{1,1}(:,1)),max(data.activity_sorted_consolidation{1,ii}(:,2))); % --> amplitude events

    if ii == 1
        data.activity_sorted_consolidation{5,ii} = zeros(size(data.activity_sorted_consolidation{2, 1},1),parameters.consolidation.total_time * parameters.sampling_rate); % --> entire session binary format
        data.activity_sorted_consolidation{6,ii} = zeros(size(data.activity_sorted_consolidation{2, 1},1),parameters.consolidation.total_time * parameters.sampling_rate); % --> entire session amplitude values
    elseif ii>1
        data.activity_sorted_consolidation{5,ii} = zeros(size(data.activity_sorted_consolidation{2, 1},1),parameters.consolidation.blocks * parameters.sampling_rate); % --> entire session binary format
        data.activity_sorted_consolidation{6,ii} = zeros(size(data.activity_sorted_consolidation{2, 1},1),parameters.consolidation.blocks * parameters.sampling_rate); % --> entire session amplitude values
    end

end

for jj = 1:size(data.events{2, 1},1)+1
    for ii = 1:size(data.activity_sorted_consolidation{1,jj},1)

        temp_data = data.activity{2,jj}(data.activity{2,jj}(:,2) == ii, 1)';

        % time in samples
        data.activity_sorted_consolidation{2,jj}(data.activity_sorted_consolidation{1, jj}(ii,1),1:length(temp_data)) = round(temp_data * parameters.sampling_rate);
        data.activity_sorted_consolidation{2,jj}(data.activity_sorted_consolidation{1, jj}(ii,1),data.activity{2,1}(ii,:) == 0) = 1; % As the data.activity was transformed from seconds to samples, times of 0 were converted to sample 1.

        % Time in seconds
        data.activity_sorted_consolidation{3,jj}(data.activity_sorted_consolidation{1, jj}(ii,1),1:length(temp_data)) = temp_data;

        % Binary events over time
        temp_data((temp_data * parameters.sampling_rate) == 0) = 1;     % As the data.activity was transformed from seconds to samples, times of 0 were converted to sample 1.
        data.activity_sorted_consolidation{5,jj}(data.activity_sorted_consolidation{1, jj}(ii,1),ceil(temp_data * parameters.sampling_rate)) = 1;  % time sample to binary

        %  Amplitude of the event
        temp_data_amplitude = data.activity{2,jj}(data.activity{2,jj}(:,2) == ii, 3)';
        data.activity_sorted_consolidation{4,jj}(data.activity_sorted_consolidation{1, jj}(ii,1),1:length(temp_data_amplitude)) = temp_data_amplitude;

        %  Amplitude of the event over time
        data.activity_sorted_consolidation{6,jj}(data.activity_sorted_consolidation{1, jj}(ii,1),ceil(temp_data * parameters.sampling_rate)) = temp_data_amplitude;  % time sample to binary

    end

    % Time vector in seconds full session
    data.activity_sorted_consolidation{7,jj} = linspace(0,size(data.activity_sorted_consolidation{6,jj},2)/parameters.sampling_rate,size(data.activity_sorted_consolidation{6,jj},2));

end

clear('temp_data','temp_data_amplitude','ii','jj')

%% Extinction - Sorting Cells activity and timestamps

% Categorizing the ID of each cell and the number of firings (events) of each one.
% GC --> numer of events
% GR --> cell ID

data.activity_sorted_extinction = [];

for ii = 1:(size(data.events{3, 1},1) + size(data.events{3, 2},1))+2
    [data.activity_sorted_extinction{1,ii}(:,2),data.activity_sorted_extinction{1,ii}(:,1)] = groupcounts(data.activity{3,ii}(:,2));
end

% In Each cell row
% -->  row: cells
% -->  columns: time

for ii = 1:(size(data.events{3, 1},1) + size(data.events{3, 2},1))+2
    data.activity_sorted_extinction{2,ii} = NaN(max(data.activity_sorted_extinction{1,1}(:,1)),max(data.activity_sorted_extinction{1,ii}(:,2))); % --> time in samples
    data.activity_sorted_extinction{3,ii} = NaN(max(data.activity_sorted_extinction{1,1}(:,1)),max(data.activity_sorted_extinction{1,ii}(:,2))); % --> time in seconds
    data.activity_sorted_extinction{4,ii} = NaN(max(data.activity_sorted_extinction{1,1}(:,1)),max(data.activity_sorted_extinction{1,ii}(:,2))); % --> amplitude events

    if ii == 1
        data.activity_sorted_extinction{5,ii} = zeros(size(data.activity_sorted_extinction{2, 1},1),parameters.extinction.total_time * parameters.sampling_rate); % --> entire session binary format
        data.activity_sorted_extinction{6,ii} = zeros(size(data.activity_sorted_extinction{2, 1},1),parameters.extinction.total_time * parameters.sampling_rate); % --> entire session amplitude values
    
    elseif ii == 2
        data.activity_sorted_extinction{5,ii} = zeros(size(data.activity_sorted_extinction{2, 1},1),parameters.extinction.baseline_time * parameters.sampling_rate); % --> entire session binary format
        data.activity_sorted_extinction{6,ii} = zeros(size(data.activity_sorted_extinction{2, 1},1),parameters.extinction.baseline_time * parameters.sampling_rate); % --> entire session amplitude values
    
    elseif ii>2 & mod(ii,2) > 0
        data.activity_sorted_extinction{5,ii} = zeros(size(data.activity_sorted_extinction{2, 1},1),parameters.extinction.CS_tone_time * parameters.sampling_rate); % --> entire session binary format
        data.activity_sorted_extinction{6,ii} = zeros(size(data.activity_sorted_extinction{2, 1},1),parameters.extinction.CS_tone_time * parameters.sampling_rate); % --> entire session amplitude values
   
    elseif ii>2 & mod(ii,2) == 0
        data.activity_sorted_extinction{5,ii} = zeros(size(data.activity_sorted_extinction{2, 1},1),parameters.extinction.ITI_time * parameters.sampling_rate); % --> entire session binary format
        data.activity_sorted_extinction{6,ii} = zeros(size(data.activity_sorted_extinction{2, 1},1),parameters.extinction.ITI_time * parameters.sampling_rate); % --> entire session amplitude values
    
    end

end

for jj = 1:(size(data.events{3, 1},1) + size(data.events{3, 2},1))+2
    for ii = 1:size(data.activity_sorted_extinction{1,jj},1)

        temp_data = data.activity{3,jj}(data.activity{3,jj}(:,2) == ii, 1)';

        % time in samples
        data.activity_sorted_extinction{2,jj}(data.activity_sorted_extinction{1, jj}(ii,1),1:length(temp_data)) = round(temp_data * parameters.sampling_rate);
        data.activity_sorted_extinction{2,jj}(data.activity_sorted_extinction{1, jj}(ii,1),data.activity{3,1}(ii,:) == 0) = 1; % As the data.activity was transformed from seconds to samples, times of 0 were converted to sample 1.

        % Time in seconds
        data.activity_sorted_extinction{3,jj}(data.activity_sorted_extinction{1, jj}(ii,1),1:length(temp_data)) = temp_data;

        % Binary events over time
        temp_data((temp_data * parameters.sampling_rate) == 0) = 1;     % As the data.activity was transformed from seconds to samples, times of 0 were converted to sample 1.
        data.activity_sorted_extinction{5,jj}(data.activity_sorted_extinction{1, jj}(ii,1),ceil(temp_data * parameters.sampling_rate)) = 1;  % time sample to binary

        %  Amplitude of the event
        temp_data_amplitude = data.activity{3,jj}(data.activity{3,jj}(:,2) == ii, 3)';
        data.activity_sorted_extinction{4,jj}(data.activity_sorted_extinction{1, jj}(ii,1),1:length(temp_data_amplitude)) = temp_data_amplitude;

        %  Amplitude of the event over time
        data.activity_sorted_extinction{6,jj}(data.activity_sorted_extinction{1, jj}(ii,1),ceil(temp_data * parameters.sampling_rate)) = temp_data_amplitude;  % time sample to binary

    end

    % Time vector in seconds full session
    data.activity_sorted_extinction{7,jj} = linspace(0,size(data.activity_sorted_extinction{6,jj},2)/parameters.sampling_rate,size(data.activity_sorted_extinction{6,jj},2));

end

clear('temp_data','temp_data_amplitude','ii','jj')

%% Retrieval - Sorting Cells activity and timestamps

% Categorizing the ID of each cell and the number of firings (events) of each one.
% GC --> numer of events
% GR --> cell ID

data.activity_sorted_retrieval = [];

for ii = 1:(size(data.events{4, 1},1) + size(data.events{4, 2},1))+2
    [data.activity_sorted_retrieval{1,ii}(:,2),data.activity_sorted_retrieval{1,ii}(:,1)] = groupcounts(data.activity{4,ii}(:,2));
end

% In Each cell row
% -->  row: cells
% -->  columns: time

for ii = 1:(size(data.events{4, 1},1) + size(data.events{4, 2},1))+2
    data.activity_sorted_retrieval{2,ii} = NaN(max(data.activity_sorted_retrieval{1,1}(:,1)),max(data.activity_sorted_retrieval{1,ii}(:,2))); % --> time in samples
    data.activity_sorted_retrieval{3,ii} = NaN(max(data.activity_sorted_retrieval{1,1}(:,1)),max(data.activity_sorted_retrieval{1,ii}(:,2))); % --> time in seconds
    data.activity_sorted_retrieval{4,ii} = NaN(max(data.activity_sorted_retrieval{1,1}(:,1)),max(data.activity_sorted_retrieval{1,ii}(:,2))); % --> amplitude events

    if ii == 1
        data.activity_sorted_retrieval{5,ii} = zeros(size(data.activity_sorted_retrieval{2, 1},1),parameters.retrieval.total_time * parameters.sampling_rate); % --> entire session binary format
        data.activity_sorted_retrieval{6,ii} = zeros(size(data.activity_sorted_retrieval{2, 1},1),parameters.retrieval.total_time * parameters.sampling_rate); % --> entire session amplitude values
    
    elseif ii == 2
        data.activity_sorted_retrieval{5,ii} = zeros(size(data.activity_sorted_retrieval{2, 1},1),parameters.retrieval.baseline_time * parameters.sampling_rate); % --> entire session binary format
        data.activity_sorted_retrieval{6,ii} = zeros(size(data.activity_sorted_retrieval{2, 1},1),parameters.retrieval.baseline_time * parameters.sampling_rate); % --> entire session amplitude values
    
    elseif ii>2 & mod(ii,2) > 0
        data.activity_sorted_retrieval{5,ii} = zeros(size(data.activity_sorted_retrieval{2, 1},1),parameters.retrieval.CS_tone_time * parameters.sampling_rate); % --> entire session binary format
        data.activity_sorted_retrieval{6,ii} = zeros(size(data.activity_sorted_retrieval{2, 1},1),parameters.retrieval.CS_tone_time * parameters.sampling_rate); % --> entire session amplitude values
   
    elseif ii>2 & mod(ii,2) == 0
        data.activity_sorted_retrieval{5,ii} = zeros(size(data.activity_sorted_retrieval{2, 1},1),parameters.retrieval.ITI_time * parameters.sampling_rate); % --> entire session binary format
        data.activity_sorted_retrieval{6,ii} = zeros(size(data.activity_sorted_retrieval{2, 1},1),parameters.retrieval.ITI_time * parameters.sampling_rate); % --> entire session amplitude values
    
    end

end

for jj = 1:(size(data.events{4, 1},1) + size(data.events{4, 2},1))+2
    for ii = 1:size(data.activity_sorted_retrieval{1,jj},1)

        temp_data = data.activity{4,jj}(data.activity{4,jj}(:,2) == ii, 1)';

        % time in samples
        data.activity_sorted_retrieval{2,jj}(data.activity_sorted_retrieval{1, jj}(ii,1),1:length(temp_data)) = round(temp_data * parameters.sampling_rate);
        data.activity_sorted_retrieval{2,jj}(data.activity_sorted_retrieval{1, jj}(ii,1),data.activity{3,1}(ii,:) == 0) = 1; % As the data.activity was transformed from seconds to samples, times of 0 were converted to sample 1.

        % Time in seconds
        data.activity_sorted_retrieval{3,jj}(data.activity_sorted_retrieval{1, jj}(ii,1),1:length(temp_data)) = temp_data;

        % Binary events over time
        temp_data((temp_data * parameters.sampling_rate) == 0) = 1;     % As the data.activity was transformed from seconds to samples, times of 0 were converted to sample 1.
        data.activity_sorted_retrieval{5,jj}(data.activity_sorted_retrieval{1, jj}(ii,1),ceil(temp_data * parameters.sampling_rate)) = 1;  % time sample to binary

        %  Amplitude of the event
        temp_data_amplitude = data.activity{4,jj}(data.activity{4,jj}(:,2) == ii, 3)';
        data.activity_sorted_retrieval{4,jj}(data.activity_sorted_retrieval{1, jj}(ii,1),1:length(temp_data_amplitude)) = temp_data_amplitude;

        %  Amplitude of the event over time
        data.activity_sorted_retrieval{6,jj}(data.activity_sorted_retrieval{1, jj}(ii,1),ceil(temp_data * parameters.sampling_rate)) = temp_data_amplitude;  % time sample to binary

    end

    % Time vector in seconds full session
    data.activity_sorted_retrieval{7,jj} = linspace(0,size(data.activity_sorted_retrieval{6,jj},2)/parameters.sampling_rate,size(data.activity_sorted_retrieval{6,jj},2));

end

clear('temp_data','temp_data_amplitude','ii','jj')


%% Plot
% 
% % Set Figure
% figure
% set(gcf,'color','white')
% box 'off'
% hold on
% 
% baseline = [1 180];
% trial    = [180 190].*parameters.sampling_rate;
% ITI      = [192 142];
% 
% for ii = 1:size(data.activity,2)
% 
%     % data.activity to plot
%     data.activity_2_plot = data.activity{5,ii};
%     time_v      = data.activity{6,ii};
% 
%     % factor
%     factor = (1:(size(data.activity{2, 1},1))).*1;
% 
%     subplot(1,2,ii)
%     % Select fields with data.activity
%     r = plot(time_v, bsxfun(@plus, data.activity_2_plot, factor'),'.','Color','[0.5, 0.5, 0.5,]','LineWidth',2);
%     hold on
%     for ii = 1:45
%         if ii ~= 1
%             r = plot(time_v(1,trial(1,1):trial(1,2)), bsxfun(@plus, data.activity_2_plot(:,trial(1,1):trial(1,2)), factor'),'.','Color','[0.6350, 0.0780, 0.1840]','LineWidth',2);
%         end
%     end
% 
%     title({'\fontsize{12}Cell Activity'; 'Extintion Session'});
%     xlabel('\fontsize{12}Time (s)');
%     ylabel('\fontsize{12}cells ');
%     xlim([0 max(time_v)]);
%     %xlim([180 190])
% end


%%
% Set Figure

% figure
% set(gcf,'color','white')
% box 'off'

baseline       = [1 data.events{2, 1}(1,1)-1];
CS_Tone_idx    = data.events{3, 1};
ITI_idx        = data.events{3, 2};

freezing = round(data.events_behavior{1, 1}./10);

for pp = 2%:size(data.activity,2)

    % data.activity to plot
    data.activity_2_plot = data.activity{5,pp};
    time_v      = data.activity{6,pp};

    subplot(1,2,pp)
    hold on

    for ii = 1:size(data.activity_2_plot,1)
        spike_ = time_v(data.activity_2_plot(ii,:)==1);
        for jj = 1:length(spike_)

            if sum(spike_(jj) >= freezing(:,1) & spike_(jj) <= freezing(:,2)) == 1
                plot([spike_(jj) spike_(jj)],[ii ii],'.','Color',[.6, 0, 0],'MarkerSize', 10)
            else
                plot([spike_(jj) spike_(jj)],[ii ii],'.','color',[0.5, 0.5, 0.5, 0.5],'MarkerSize', 10)

            end

        end
        
        plot([time_v(data.events{2, 1}(:,1));time_v(data.events{2, 1}(:,2))], 200.*([ones(1,length(data.events{2, 1}(:,1))).*.5;ones(1,length(data.events{2, 1}(:,2))).*.5]),'k-','linew', 20,'Color',[1, .4, .4])
        %ylim([0 100])


        plot(time_v, movmean(data.behavior{1, 1}(1:22431),40)./10-10,'linew', 1, 'Color',[0.6, 0.6, 0.6, 0.4])

        ylim([-10 100])

    end



    xlim([0 max(time_v)]);
    %ylim([0 100]);

    xlabel('Time (s)')
    ylabel('#Cells (1-90)')
    a = gca; % Get axis
    %a.YColor = 'k';
    a.YTick = []; 

title('Retrieval')

end








