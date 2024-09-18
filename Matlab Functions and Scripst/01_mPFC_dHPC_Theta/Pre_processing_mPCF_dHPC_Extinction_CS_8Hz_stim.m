
%% Data pre-processing

% (First extract the data using the "Extracting raw LFPs and Events" script)

% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  12/2023
% Last update: 05/2024

%% First extract the data with : Extracting_LFPs_and_events.m

%%
fprintf('\n Data pre-processing... \n');
%%
% - CHANNELS MAP - Monopolar -

%  Design
%  Ch 1-8 was always mPFC IL
%  Ch 9-16 was always mPFC PL
%  Ch 17-32 was always dHPC

%% The data collected has already been downsampled to 1kHz.

% Data downsampling
% parameters.downsampling = 40;
% 
% for ii = 1:size(data.lfp{1,1},1)
%     data.lfp{5,1}(ii,:) = decimate(data.lfp{1,1}(ii,:),parameters.downsampling);
% end

parameters.decimated_srate = parameters.original_srate; %./parameters.downsampling;
data.timev_decimated = linspace(0,length(data.lfp{1,1})/parameters.decimated_srate,length(data.lfp{1,1}));

%% The data collected has already been downsampled to 1kHz.

% Here, I will only reorganize the variable Data.lfp.
% In other words, downsampled data will be placed on line 5. The first cells will eventually be designated for raw records, meaning with the original sampling frequency.
% Furthermore, in this session, the channels for analysis will be selected.


% Time Session for Contextual Exposure. Normalizing data to the same time window
% session_length = 1200; %s
% 
% if length(data.lfp{1,1}) > session_length*parameters.original_srate
% 
%     data.lfp{5,1} = data.lfp{1,1}([5,13,24],1:(session_length*parameters.original_srate));
%     data.lfp{1,1} = [];
% else
%     sample2add = session_length*parameters.original_srate - length(data.lfp{1,1});
%     data.lfp{5,1} = [data.lfp{1,1}([5,13,24],:) zeros(size(data.lfp{1,1}([5,13,24],:),1),sample2add)];
%     data.lfp{1,1} = [];
% 
% end


% Time Session
data.lfp{5,1} = data.lfp{1,1};
data.lfp{1,1} = [];

%% Use first channel from each probe as Reference
% temp1 = data.lfp{5, 1}(1:8,:) - data.lfp{5, 1}(1,:);    % mPFC IL
% temp2 = data.lfp{5, 1}(9:16,:) - data.lfp{5, 1}(11,:);   % mPFC PL
% temp3 = data.lfp{5, 1}(17:32,:) - data.lfp{5, 1}(18,:); % dHPC
% 
% data.lfp{5, 1} = [];
% data.lfp{5, 1} = [temp1;temp2;temp3];
% 
% clear('temp1','temp2')

%% Filter Data
% Set of filters by VRCarva (https://github.com/vrcarva) and EEGLab

% % Define frequencies cutoff
% 
% % fun_myfilters function - NOTCH
% parameters.filter.notch1          = [58 62];
% parameters.filter.notch2          = [118 122];     
% parameters.filter.notch3          = [178 182];
% 
% % parameters.filter.longcutoff_1      = [1 0];      %2
% % parameters.filter.longcutoff_2      = [0 300];    %2
% % parameters.filter.thetacutoff_1     = [2 12];     %3
% % parameters.filter.thetacutoff_2     = [3 6];      %4
% % parameters.filter.thetacutoff_3     = [6 9];      %5
% % parameters.filter.lowgammacutoff    = [30 50];    %6
% % parameters.filter.highgammacutoff   = [62 100];   %7
% 
% % Notch Filter - on(1) /off(0)
% params.bandstop = 1;
% 
% % Each cell --> columns: parameters.filters according to the above order 
% 
% % Filter decimate data
% for jj = 1:size(data.lfp{5,1},1)
% 
%     data.data{5,2}(jj,:)  = fun_myfilters(data.lfp{5,1}(jj,:),parameters.decimated_srate,parameters.filter.notch1,'iir',params);
%     data.data{5,2}(jj,:)  = fun_myfilters(data.lfp{5,2}(jj,:),parameters.decimated_srate,parameters.filter.notch2 ,'iir',params);
%     data.data{5,2}(jj,:)  = fun_myfilters(data.lfp{5,2}(jj,:),parameters.decimated_srate,parameters.filter.notch3,'iir',params);
% 
% %   data.lfp{5,2}(jj,:)  = fun_myfilters(data.lfp{5,1}(jj,:),parameters.decimated_srate,parameters.filter.longcutoff_1,'eegfilt',params);
% %   data.lfp{5,2}(jj,:)  = fun_myfilters(data.lfp{5,2}(jj,:),parameters.decimated_srate,parameters.filter.longcutoff_2,'eegfilt',params);
% %     
% %   data.lfp{5,3}(jj,:)  = fun_myfilters(data.lfp{5,1}(jj,:),parameters.decimated_srate,parameters.filter.thetacutoff_1,'eegfilt',params);
% %   data.lfp{5,4}(jj,:)  = fun_myfilters(data.lfp{5,1}(jj,:),parameters.decimated_srate,parameters.filter.thetacutoff_2,'eegfilt',params);
% %   data.lfp{5,5}(jj,:)  = fun_myfilters(data.lfp{5,1}(jj,:),parameters.decimated_srate,parameters.filter.thetacutoff_3,'eegfilt',params);
% %
% %   data.lfp{5,6}(jj,:)  = fun_myfilters(data.lfp{5,1}(jj,:),parameters.decimated_srate,parameters.filter.lowgammacutoff,'eegfilt',params);
% %   data.lfp{5,7}(jj,:)  = fun_myfilters(data.lfp{5,1}(jj,:),parameters.decimated_srate,parameters.filter.highgammacutoff,'eegfilt',params);
% % 
% end
% 
% clear('jj','params')

%% Set CS and ITI timestamps
% data.events
% cell line 1 : raw events from Plexon

% cell line 2 : STM and CS-Tones
%    __ columns 1 --> 1: STM and CS-Tones
%    __ columns 2 --> 1: ITI

% cell line 3 : CS-Tones
%    __ columns 1 --> 1: STM and CS-Tones



% % ___________________________________________________________
% % ---- IMPORTANT ----
% %
% % The timestamps were decimated along with the records.
% % ___________________________________________________________


% Extract timestamps
dataascell = struct2cell(data.events{1, 1});
emptyCells = cellfun(@isempty,dataascell(2,:));

% Case 1: it will be add 10 exact seconds from the timestamp to define the stimulus window 
% and  30 seconds backward from the timestamps to define the ITI window. 
% Note that a small time window will be add to the CS-trial and a small window will be lost to the ITI. 

% Time in samples
data.events{2,1}(:,1) = cell2mat(dataascell(1,emptyCells==0))'; % STM and CS-Start

% Separate each event
event_temp1 = data.events{2,1}(1:2:end,1); % STM Start
event_temp2 = data.events{2,1}(2:2:end,1); % CS-Tone Start
data.events{2,1} = [];
data.events{2,1} = event_temp1;
data.events{3,1} = event_temp2;

% Define Trials STM + CS-Tones
ST_CS_time               = 20*parameters.original_srate;           % Defining ST_CS time
data.events{2,1}(:,2)    = data.events{2,1}(:,1)+ST_CS_time;       % ST_CS-end

% Define Trials STM + CS-Tones
CS_time                  = 10*parameters.original_srate;           % Defining CS time
data.events{3,1}(:,2)    = data.events{3,1}(:,1)+CS_time ;         % CS-end

% Define ITI-Trials
ITI_time              = 30*parameters.original_srate;               % Defining ITI time

data.events{2,2}(:,1) = data.events{2,1}(2:end,1) - ITI_time -1;    % ITI-Start
data.events{2,2}(:,2) = data.events{2,1}(2:end,1) - 1;              % ITI-Stop

....% last ITI-Trial. Normalize to 30 s after CS-Trial
data.events{2,2}(end+1,1) = data.events{2,1}(end,2) +1 ;            % ITI-Start
data.events{2,2}(end,2)   = data.events{2,1}(end,2) + ITI_time + 1; % ITI-Stop


% Plot to confirm
% figure
% plot([data.timev_decimated(data.events{2, 1}(:,1));data.timev_decimated(data.events{2, 1}(:,2))], [ones(1,length(data.events{2, 1}(:,1)));ones(1,length(data.events{2, 1}(:,2)))],'k-','linew', 2,'Color',[1, .4, .4])
% ylim([.5 1.5])


% CS and ITI - Time in sec.
% data.events{3,1} = data.events{2,1}./parameters.original_srate;
% data.events{3,2} = data.events{2,2}./parameters.original_srate;


clear('locs','TS_temp','TS_temp_1','TS_temp_2','v','ITI_time','CS_time')

%% Organizing data - Considering only the baseline CS and ITI Events

% Cell rows
% 1 - 4: empty (future units Analysis)
% 5 : Decimate data
% 6 : baseline
% 7 : CS-TONES
% 8 : ITI

% Cell columns --> frequencies cutoff according filter data
% in each cell
% - rows        - > channels
% - columns     - > time
% - 3 dimension - > behavioral events; CS-TONES or ITIs.


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

% STM and CS-Tones
for ii = 1:size(data.lfp,2)
    for jj = 1:size(data.events{2, 1},1)

        data.lfp{7,ii}(:,:,jj) = data.lfp{5,ii}(:,(round((data.events{2, 1}(jj,1))./parameters.downsampling)):(round((data.events{2, 1}(jj,2))./parameters.downsampling)-1));

        if isempty(data.lfp{1,ii}) % Future raw data - The raw data has never been filtered until now.
            continue
        else
%             data.lfp{3,ii}(:,:,jj) = data.lfp{1,ii}(:,data.events{3, 1}(jj,1):data.events{3, 1}(jj,2)-1);
        end
    end
end

% ITI
for ii = 1:size(data.lfp,2)
    for jj = 1:size(data.events{2, 2},1)
            
        data.lfp{8,ii}(:,:,jj) = data.lfp{5,ii}(:,(round((data.events{2, 2}(jj,1))./parameters.downsampling)):(round((data.events{2, 2}(jj,2))./parameters.downsampling)-1));

        if isempty(data.lfp{1,ii}) % Future raw data - The raw data has never been filtered until now.
            continue
        else
%             data.lfp{4,ii}(:,:,jj) = data.lfp{1,ii}(:,data.events{3, 2}(jj,1):data.events{3, 2}(jj,2)-1);
        end
    end
end


% CS-Tones
for ii = 1:size(data.lfp,2)
    for jj = 1:size(data.events{3, 1},1)

        data.lfp{9,ii}(:,:,jj) = data.lfp{5,ii}(:,(round((data.events{3, 1}(jj,1))./parameters.downsampling)):(round((data.events{3, 1}(jj,2))./parameters.downsampling)-1));

        if isempty(data.lfp{1,ii}) % Future raw data - The raw data has never been filtered until now.
            continue
        else
%             data.lfp{3,ii}(:,:,jj) = data.lfp{1,ii}(:,data.events{4, 1}(jj,1):data.events{4, 1}(jj,2)-1);
        end
    end
end


clear('ii','jj')

%% Behavior analysis

Behavior_mPCF_dHPC_Extinction_CS_8Hz_stim

%% Set index of freezing and non-freezing for baseline, CS-Tones or ITI events

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


% Initialize
data.events_behavior =[];


% Baseline
% Freezing idx
baseline_idx_f = data.behavior{2,1}(1,:) < data.events{2, 1}(1,1);  % logical idx
data.events_behavior{1,1}(:,1) = data.behavior{2,1}(1,baseline_idx_f);                                           % Start idx
data.events_behavior{1,1}(:,2) = data.behavior{2,1}(1,baseline_idx_f)'  + data.behavior{2,1}(2,baseline_idx_f)';  % End idx

% Non-freezing
baseline_idx_nf = data.behavior{2,2}(1,:) < data.events{2, 1}(1,1); % logical idx
data.events_behavior{2,1}(:,1) = data.behavior{2,2}(1,baseline_idx_nf);                                           % Start idx
data.events_behavior{2,1}(:,2) = data.behavior{2,2}(1,baseline_idx_nf)' + data.behavior{2,2}(2,baseline_idx_nf)' - 1; % End idx


% CS and ITI
for ii = 1:size(data.events{2, 1},1)
    
    % CS-Tone
    % Freezing idx
    CS_idx = data.behavior{2,1}(1,:) >= data.events{2, 1}(ii,1) & data.behavior{2,1}(1,:) <= data.events{2, 1}(ii,2);   % logical idx inside CS-TONE
    data.events_behavior{3,ii}(:,1) = data.behavior{2,1}(1,CS_idx);                                       % Start idx
    data.events_behavior{3,ii}(:,2) = data.behavior{2,1}(1,CS_idx) + data.behavior{2,1}(2,CS_idx) - 1;    % End idx
    % Non-freezing
    CS_idx_1 = data.behavior{2,2}(1,:) >= data.events{2, 1}(ii,1) & data.behavior{2,2}(1,:) <= data.events{2, 1}(ii,2); % logical idx inside CS-TONE
    data.events_behavior{4,ii}(:,1) = data.behavior{2,2}(1,CS_idx_1);                                  % Start idx
    data.events_behavior{4,ii}(:,2) = data.behavior{2,2}(1,CS_idx_1) + data.behavior{2,2}(2,CS_idx_1) - 1; % End idx

    % ITI
    % Freezing idx    
    ITI_idx = data.behavior{2,1}(1,:) >= data.events{2, 2}(ii,1) & data.behavior{2,1}(1,:) <= data.events{2, 2}(ii,2); % logical idx inside ITI
    data.events_behavior{5,ii}(:,1) = data.behavior{2,1}(1,ITI_idx);                                 % Start idx
    data.events_behavior{5,ii}(:,2) = data.behavior{2,1}(1,ITI_idx) + data.behavior{2,1}(2,ITI_idx) - 1; % End idx
    % Non-freezing
    ITI_idx = data.behavior{2,2}(1,:) >= data.events{2, 2}(ii,1) & data.behavior{2,2}(1,:) <= data.events{2, 2}(ii,2); % logical idx inside ITI
    data.events_behavior{6,ii}(:,1) = data.behavior{2,2}(1,ITI_idx);                                 % Start idx
    data.events_behavior{6,ii}(:,2) = data.behavior{2,2}(1,ITI_idx) + data.behavior{2,2}(2,ITI_idx) - 1; % End idx

end

clear('baseline_idx_f','baseline_idx_nf','CS_idx','CS_idx_1','ITI_idx','ii')

%% Organizing data - Considering freezing and non-freezing behavior


% Plot to check the best time win length
% From the histogram, select the time window that encompasses the largest number of events.
% h_freezing    = ((data.behavior{2,1}(1,:)+data.behavior{2,1}(2,:)) - data.behavior{2,1}(1,:))./parameters.original_srate;
% h_Nonfreezing = ((data.behavior{2,2}(1,:)+data.behavior{2,2}(2,:)) - data.behavior{2,2}(1,:))./parameters.original_srate;
% 
% histogram(h_freezing,100)
% histogram(h_Nonfreezing,100)



% First Option
% Select events considering a maximum time windows
% timewin = 5; % sec.
% 
% for ii = 1:size(data.events_behavior,1)
%     for jj = 1:size(data.events_behavior,2) 
% 
%         if isempty(data.events_behavior{ii, jj})
%             continue
%         end
% 
%         f_idx{ii,jj} = data.events_behavior{ii, jj}(((data.events_behavior{ii, jj}(:,2) - data.events_behavior{ii, jj}(:,1))<=timewin.*parameters.decimated_srate),:);
% 
%     end
% end



% Second Option
% Define time win for zero pad or to cut the data

timewin = 5; % sec.
f_idx = data.events_behavior; %Redundant, but it makes the code cleaner


% data.events_behavior organization - From data.lfp{5,1} - full decimated data
% - Row 1: Baseline freezing
% - Row 2: Baseline not freezing
% - Row 3: CS-TONE freezing
% - Row 4: CS-TONE not freezing
% - Row 5: ITI freezing
% - Row 6: ITI not freezing

%   - Columns trials

%       Within each cell:
%           - Rows            : channels
%           - Columns         : time
%           - third dimension : events

% ---------------------
% IMPORTANTE NOTE:
% The data was zero-padded to occupy a maximum size of timewin seconds, so normalizing the frequency resolution for the spectra. 
% Events lasting more than timewin seconds were set to exactly 5 seconds.
% ---------------------



data.lfp_behavior = [];

for ii = 1:size(f_idx,1,1)
    for jj = 1:size(f_idx,2)
        
        if isempty(f_idx{ii,jj})
            continue
        end
        
        data.lfp_behavior{ii,jj} = zeros(size(data.lfp{5,1},1), timewin * parameters.decimated_srate, size(f_idx{ii,jj},1)); % fixing timmig win...zero padding

        for ff = 1:size(f_idx{ii,jj},1)

            if (f_idx{ii,jj}(ff,2) - f_idx{ii,jj}(ff,1) >= timewin * parameters.decimated_srate) % Conditional for Time events >= than timewin
               data.lfp_behavior{ii,jj}(:,1:size(data.lfp_behavior{ii,jj},2),ff) = data.lfp{5,1}(:,f_idx{ii,jj}(ff,1):(f_idx{ii,jj}(ff,1) + timewin * parameters.decimated_srate)-1);
            else
               
               temp_data = data.lfp{5,1}(:,f_idx{ii,jj}(ff,1):f_idx{ii,jj}(ff,2)); % events < time_win
               data.lfp_behavior{ii,jj}(:,1:size(temp_data,2),ff) = temp_data;

            end
        end            
    end    
end

clear('ii','jj','ff','timewin','f_idx','m_idx','time_win_pos','time_win_pos','time_win_pre','temp_data' )

%% Plot baseline to remove noise epochs
% Checking through the hippocampal channel, which generally exhibits more noise.

%freezing
% figure
% for ii = 1:size(data.lfp_behavior{1,1},3)
%     subplot(3,3,ii)
%     plot(data.lfp_behavior{1,1}(3,:,ii))
% end
% 
% %Non-freezing
% figure
% for ii = 1:size(data.lfp_behavior{2,1},3)
%     subplot(1,1,ii)
%     plot(data.lfp_behavior{2,1}(3,:,ii))
% end
% 
% clear('ii')

% Exclude baseline noise events
%data.lfp_behavior{1,1}(:,:,[2]) = [];
%data.lfp_behavior{2,1}(:,:,[1]) = [];

%% Save

% newStr = regexprep(files.id.name,'.mat','_');
% newStr = files.id(ms).name(1:end-8);
newStr = id(1:end-5);

path = '/Users/flavio/Desktop';
% path = files.FilesLoaded{1, 1}.folder;
% 
name = strcat(path,'/',newStr);

save(name,'data','parameters','id','-v7.3')

clear('newStr','path')

%% last update 19/05/2024 
%  listening: 
