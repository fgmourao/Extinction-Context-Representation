
%% Data pre-processing

% (First extract the data using the "Extracting raw LFPs and Events" script)

% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  12/2023
% Last update: 01/2024

%%
% - CHANNELS MAP - Monopolar -
% - CHANNELS MAP - Bipolar -

%% First extract the data with : Extracting_LFPs_and_events.m

%%

fprintf('\n Data pre-processing... \n');


%% Delete reference channel

data.lfp{1,1}(1,:) = [];

%% Remove bad channels

% Animal 0 - habituation
data.lfp{1,1}([1 4 5 8 12],:,:) = [];

%% Bipolar derivation

data.lfp{1,1} = diff(data.lfp{1,1},1,1);

%% Data downsampling

parameters.downsampling = 40;

for ii = 1:size(data.lfp{1,1},1)
    data.lfp{5,1}(ii,:) = decimate(data.lfp{1,1}(ii,:),parameters.downsampling);
end

parameters.decimated_srate = parameters.original_srate./parameters.downsampling;
data.timev_decimated = linspace(0,length(data.lfp{5,1})/parameters.decimated_srate,length(data.lfp{5,1}));

clear('ii')

%% Filter Data
% Set of filters by VRCarva (https://github.com/vrcarva) and EEGLab

% Define frequencies cutoff

% fun_myfilters function - NOTCH
% parameters.filter.notch1          = [58 62];
% parameters.filter.notch2          = [118 122];     
% parameters.filter.notch3          = [178 182];

parameters.filter.longcutoff_1        = [1 0];    %2
parameters.filter.longcutoff_2        = [0 140];    %2

parameters.filter.thetacutoff_1     = [2 12];     %3

% Notch Filter - on(1) /off(0)
params.bandstop = 0;

% Each cell --> columns: parameters.filters according to the above order 

% Filter decimate data
for jj = 1:size(data.lfp{5,1},1)

    data.lfp{5,2}(jj,:)  = fun_myfilters(data.lfp{5,1}(jj,:),parameters.decimated_srate,parameters.filter.longcutoff_1,'eegfilt',params);
    data.lfp{5,2}(jj,:)  = fun_myfilters(data.lfp{5,1}(jj,:),parameters.decimated_srate,parameters.filter.longcutoff_2,'eegfilt',params);

    data.lfp{5,3}(jj,:)  = fun_myfilters(data.lfp{5,1}(jj,:),parameters.decimated_srate,parameters.filter.thetacutoff_1,'eegfilt',params);

end

clear('jj','params')


%% Set CS and ITI timestamps

% Tugce experiment. 
% Current configuration. CS is a sine wave. Timestamps only reference the:
% Start experiment. 
% CS - start.
% End experiment

% Considering that we only have timestamps related to the start periods e stop session, it will be add 10 exact seconds from the timestamp to define the stimulus window 
% and  30 seconds backward from the timestamps to define the ITI window. 
% Note that a small time window will be lost. However, this will ensure precision regarding the event.



% IMPORTANT: 
% In the record made by Tuge, in the current record setup, some timestamps may be missing. It is necessary to check and correct manually.
TS_temp = [data.events{1, 1}.sample];
TS_temp_1 = diff(TS_temp)./parameters.original_srate
TS_temp_2 = TS_temp./parameters.original_srate



% Return after proper corrections
TS_temp = TS_temp';

% Define CS-Trials

CS_time = 10*parameters.original_srate; % Defining CS time

% Time in samples
data.events{2,1}(:,1) = TS_temp(1:end);                     % CS-Start
data.events{2,1}(:,2) = TS_temp(1:end)+CS_time;             % CS-end

% Define ITI-Trials

% IMPORTANT.
% Last line in data.events correspond to the STOP timestamp. There is, when
% the experiment ends


% Case to calculate the precise timing of ITI.
%ITI_time = (data.events{2,1}(2:end,1) - data.events{2,1}(1:end-1,2))-2; % Defining ITI time


ITI_time           = 30*parameters.original_srate; % Defining ITI time
ITI_time_lasttrial = 60*parameters.original_srate; % Defining ITI time for the last trial

% Samples
data.events{2,2}(:,1) = data.events{2,1}(2:end-1,1) - ITI_time -1;       % ITI-Start
data.events{2,2}(:,2) = data.events{2,1}(2:end-1,1) - 1;              % ITI-Stop

....% last ITI-Trial. Normalize to 30 s after CS-Trial
data.events{2,2}(end+1,1) = data.events{2,1}(end,1) - ITI_time_lasttrial ;       % ITI-Start
data.events{2,2}(end,2)   = data.events{2,2}(end,1) + ITI_time;                    % ITI-Stop


% CS and ITI - Time in sec.
data.events{3,1} = data.events{2,1}./parameters.original_srate;
data.events{3,2} = data.events{2,2}./parameters.original_srate;


% Remove last line to correct. 
% A bit confusing, but it has been kept until here only for the precise calculation of the ITI interval so far.
data.events{2,1}(end,:) = []; 
data.events{3,1}(end,:) = []; 


clear('CS_time','ii','ITI_time','ITI_time_lasttrial','TS_temp','TS_temp_1','TS_temp_2')

%% Organizing data - Considering only the Events

% Cell columns --> frequencies cutoff according filter data
% in each cell
% - rows        - > channels
% - columns     - > time
% - 3 dimension - > behavioral events

% Cutting time windows

% Baseline
for ii = 1:size(data.lfp,2)

    data.lfp{6,ii} = data.lfp{5,ii}(:,1:(round((data.events{2, 1}(1,1)-1)./parameters.downsampling)));

    if isempty(data.lfp{1,ii}) % The raw data has never been filtered until now.
        continue
    else
        data.lfp{2,ii} = data.lfp{1,ii}(:,1:data.events{2, 1}(1,1)-1);
    end

end

% CS-Trials
for ii = 1:size(data.lfp,2)
    for jj = 1:size(data.events{2, 1},1)

        data.lfp{7,ii}(:,:,jj) = data.lfp{5,ii}(:,(round((data.events{2, 1}(jj,1))./parameters.downsampling)):(round((data.events{2, 1}(jj,2))./parameters.downsampling)-1));

        if isempty(data.lfp{1,ii}) % The raw data has never been filtered until now.
            continue
        else
            data.lfp{3,ii}(:,:,jj) = data.lfp{1,ii}(:,data.events{2, 1}(jj,1):data.events{2, 1}(jj,2)-1);
        end
    end
end

% ITI-Trials
for ii = 1:size(data.lfp,2)
    for jj = 1:size(data.events{2, 1},1)

        data.lfp{8,ii}(:,:,jj) = data.lfp{5,ii}(:,(round((data.events{2, 2}(jj,1))./parameters.downsampling)):(round((data.events{2, 2}(jj,2))./parameters.downsampling)-1));

        if isempty(data.lfp{1,ii}) % The raw data has never been filtered until now.
            continue
        else
            data.lfp{4,ii}(:,:,jj) = data.lfp{1,ii}(:,data.events{2, 2}(jj,1):data.events{2, 2}(jj,2)-1);
        end
    end
end

%% Run behavior analysis Script
% Behavior

%% last update 10/01/2024 
%  listening: 
