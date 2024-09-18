%% Analysis of cellular activity RAW SIGNAL. Calcium imaging experiment. Miniscope.
%  Pre-Processing - Sorting Session, events and cells

%  After preprocessing, the data.raw from both sessions_activity were merged and saved in the same spreadsheet.
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
data_excel_raw = Rat1cellsraw; % Spreedsheet name

% Let`s sorting by time to find the gap between sessions_activity and thus separate them
time_sort_raw = sortrows(data_excel_raw,1,'ascend');

% index where each session started. 
% Considering that each session lasts a minimum of 300 seconds, 250 was a reasonable number.
session_idx_ = find([1 ;diff(time_sort_raw.VarName1) > 250]);

data.raw = {};

for jj = 1:length(session_idx_)
    for ii = 1:size(time_sort_raw,2)

        % Session

        if jj == length(session_idx_)
            data.raw{jj,1}(:,ii) = table2array(time_sort_raw(session_idx_(jj):end,ii));

        else
            data.raw{jj,1}(:,ii) = table2array(time_sort_raw(session_idx_(jj):session_idx_(jj+1)-1,ii));

        end
    end
end


% Correcting the time for each session
for ii = 1:size(data.raw,1)
    data.raw{ii,1}(:,1) = data.raw{ii,1}(:,1) - data.raw{ii,1}(1,1);
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

% data.raw

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


% Baseline
for ii = 1:size(data.events{1, 1},1)
    data.raw{1,ii+1} = data.raw{1,1}(data.raw{1,1}(:,1) >= data.events{1, 1}(ii,1) & data.raw{1,1}(:,1) <= data.events{1, 1}(ii,2),:);
    
    if ii > 1
        data.raw{1,ii+1}(:,1) = data.raw{1,ii+1}(:,1) - data.events{1, 1}(ii,1);
    end
    
end


% Consolidation
for ii = 1:size(data.events{2, 1},1)
    data.raw{2,ii+1} = data.raw{2,1}(data.raw{2,1}(:,1) >= data.events{2, 1}(ii,1) & data.raw{2,1}(:,1) <= data.events{2, 1}(ii,2),:);
    
    if ii > 1
        data.raw{2,ii+1}(:,1) = data.raw{2,ii+1}(:,1) - data.events{2, 1}(ii,1);
    end

end


% Extinction
% Baseline
data.raw{3,2} = data.raw{3,1}(data.raw{3,1}(:,1) < parameters.extinction.baseline_time,:);
% CS-TONES
for ii = (1:2:(size(data.events{3, 1},1) + size(data.events{3, 2},1)))+2
    data.raw{3,ii} = data.raw{3,1} (data.raw{3,1}(:,1) >= data.events{3, 1}(floor(ii/2),1) & data.raw{3,1}(:,1) <= data.events{3, 1}(floor(ii/2),2),:);
    data.raw{3,ii+1} = data.raw{3,1}(data.raw{3,1}(:,1) >= data.events{3, 2}(floor(ii/2),1) & data.raw{3,1}(:,1) <= data.events{3, 2}(floor(ii/2),2),:);
end

% Correcting time for each trial
% CS-Tone
for ii = 3:2:(parameters.extinction.number_CS_tones+parameters.extinction.number_CS_tones)+2
        data.raw{3,ii}(:,1) = data.raw{3,ii}(:,1) - data.events{3, 1}(floor(ii/2),1);
end

% ITI
for ii = 4:2:(parameters.extinction.number_ITI+parameters.extinction.number_ITI)+2
        data.raw{3,ii}(:,1) = data.raw{3,ii}(:,1) - data.events{3, 2}((ii/2)-1,1);
end


% Retrieval
% Baseline
data.raw{4,2} = data.raw{4,1}(data.raw{4,1}(:,1) < parameters.retrieval.baseline_time,:);
% CS-TONES
for ii = (1:2:(size(data.events{4, 1},1) + size(data.events{4, 2},1)))+2
    data.raw{4,ii} = data.raw{4,1} (data.raw{4,1}(:,1) >= data.events{4, 1}(floor(ii/2),1) & data.raw{4,1}(:,1) <= data.events{4, 1}(floor(ii/2),2),:);
    data.raw{4,ii+1} = data.raw{4,1}(data.raw{4,1}(:,1) >= data.events{4, 2}(floor(ii/2),1) & data.raw{4,1}(:,1) <= data.events{4, 2}(floor(ii/2),2),:);

end

% Correcting time for each trial
% CS-Tone
for ii = 3:2:(parameters.extinction.number_CS_tones+parameters.extinction.number_CS_tones)+2
        data.raw{4,ii}(:,1) = data.raw{4,ii}(:,1) - data.events{4, 1}(floor(ii/2),1);
end

% ITI
for ii = 4:2:(parameters.extinction.number_ITI+parameters.extinction.number_ITI)+2
        data.raw{4,ii}(:,1) = data.raw{4,ii}(:,1) - data.events{4, 2}((ii/2)-1,1);
end

clear('ii')

%% Plots
data2plot = data.raw{1,1}(:,2:end)';
%data2plot(:,isnan(data2plot)) = [];
[data2plot_normalized, PS]=mapminmax(data2plot,0,1);
timev = data.raw{1,1}(:,1)';

% norm_data = (bla - minVal) / ( maxVal - minVal )
% your_original_data = minVal + norm_data.*(maxVal - minVal)

% x_normalized_z = zscore(x_normalized,0,1);

%%

figure
set(gcf,'color','white')
box 'off'
hold on
% Choose channels to plot
channels = 1:57;

% factor
factor = (channels)'*1;
r = plot(timev, bsxfun(@plus, data2plot_normalized(channels,:), factor),'Color','[0.3, 0.3, 0.3]','linew',1);
a = gca;
a.YColor = 'w';
a.YTick = [];
a.XLim = [0  timev(end)];
box off
%xline(180)
%xlim([0 300])

%%
figure
colormap parula
contourf(timev,1:57,data2plot_normalized(1:57,:),80,'linecolor','none');
colorbar
xlim([150 165])
%clim([0 1])


