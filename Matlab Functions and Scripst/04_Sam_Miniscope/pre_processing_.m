%% Load data.activity
% After preprocessing, the data.activity from both extinction and retrieval sessions were mixed and saved in the same spreadsheet.
% the table includes time (seconds) where the event occurs, the cell ID that
% produced the event, and the amplitude (value) of each event.

% The MedPC trigger starts the recording, so 0s on the miniscope recording is also 0s on the freezing output


%% Load and sorting Extinction and Retrieval sessions

% Let`s sorting by time to find the gap between sessions and thus separate them
time_sort = sortrows(SSC6preprocesseddata,1,'ascend');
[~ , idx_] = max(diff(time_sort.Times)); % index where extincton end and retrieval start

% column 1 --> extinction
% Column 2 --> retrieval

data.activity = {};

for ii = 1:size(time_sort,2)

    % Extinction
    data.activity{1,1}(:,ii) = table2array(time_sort(1:idx_,ii));
    % Retrieval
    data.activity{1,2}(:,ii) = table2array(time_sort(idx_+1:end,ii));

end

% Correcting the time of retrieval
data.activity{1,2}(:,1) = data.activity{1,2}(:,1) - data.activity{1,2}(1,1);


clear('idx_','ii','SSC6preprocesseddata','time_sort')

%% Settings

parameters.sampling_rate = 10;

parameters.baseline_time  = 180 * parameters.sampling_rate; %samples
parameters.CS_tone_time   = 10  * parameters.sampling_rate; %samples
parameters.ITI_time       = 30  * parameters.sampling_rate; %samples



%% data.events


data.events = [];
data.events{1,1} = (1:parameters.CS_tone_time + parameters.ITI_time:(parameters.CS_tone_time + parameters.ITI_time)*45) + parameters.baseline_time ;

% CS_Tones
data.events{2,1}(:,1) = data.events{1,1}';                                     % Begin in samples
data.events{2,1}(:,2) = (data.events{2,1}(:,1) + parameters.CS_tone_time) - 1; % End in samples
data.events{3,1} = data.events{2,1}./parameters.sampling_rate;                 % Begin/end in seconds

% ITI
data.events{2,2}(:,1) = (data.events{1,1}(2:end)' - parameters.ITI_time); % Begin
data.events{2,2}(:,2) = data.events{1,1}(2:end)' - 1;                     % End

data.events{2,2}(end+1,1) = data.events{2, 1}(end,2) + 1                  % Begin the last ITI epoch
data.events{2,2}(end,2) = data.events{2, 1}(end,2) + parameters.ITI_time; % End the last ITI epoch

data.events{3,2} = data.events{2,2}./parameters.sampling_rate;            % Begin/End in seconds


%% Sorting Cells and timestamps

% Extinction

% Categorizing the ID of each cell and the number of firings (events) of each one.
% GC --> numer of events
% GR --> cell ID

% Important, some cells show no activity in the extinction session, only in retrieval. And vice versa. Therefore, some IDs may not appear.

[GC_E,GR_E] = groupcounts(data.activity{1,1}(:,2));

% In Each cell row
% -->  row: cells
% -->  columns: time

data.activity{2,1} = NaN(max(GR_E),max(GC_E)); % --> time in samples
data.activity{3,1} = NaN(max(GR_E),max(GC_E)); % --> time in seconds
data.activity{4,1} = NaN(max(GR_E),max(GC_E)); % --> amplitude events

data.activity{5,1} = zeros(size(data.activity{2, 1},1),round(data.activity{1, 1}(end,1).*parameters.sampling_rate)); % --> entire session binary format
data.activity{6,1} = zeros(size(data.activity{2, 1},1),round(data.activity{1, 1}(end,1).*parameters.sampling_rate)); % --> entire session amplitude values


for ii = 1:max(data.activity{1,1}(:,2))

    % Time in samples
    temp_data.activity = data.activity{1,1}(data.activity{1,1}(:,2) == ii, 1)';
    data.activity{2,1}(ii,1:length(temp_data.activity)) = round(temp_data.activity * parameters.sampling_rate); % time in samples
    data.activity{2,1}(ii,data.activity{2,1}(ii,:) == 0) = 1; % As the data.activity was transformed from seconds to samples, times of 0 were converted to sample 1.

    % Time in seconds
    data.activity{3,1}(ii,1:length(temp_data.activity)) = temp_data.activity;

    % Binary events over time
    temp_data.activity((temp_data.activity * parameters.sampling_rate) == 0) = 1;     % As the data.activity was transformed from seconds to samples, times of 0 were converted to sample 1.
    data.activity{5,1}(ii,round(temp_data.activity * parameters.sampling_rate)) = 1;  % time sample to binary

    %  Amplitude of the event
    temp_data.activity_amplitude = data.activity{1,1}(data.activity{1,1}(:,2) == ii, 3)';
    data.activity{4,1}(ii,1:length(temp_data.activity_amplitude)) = temp_data.activity_amplitude;

    %  Amplitude of the event over time
    data.activity{6,1}(ii,round(temp_data.activity * parameters.sampling_rate)) = temp_data.activity_amplitude;  % time sample to binary

end

% Time vector in seconds full session
data.activity{7,1} = linspace(0,data.activity{1, 1}(end,1),round(data.activity{1, 1}(end,1).*parameters.sampling_rate));



% Retrieval

% Categorizing the ID of each cell and the number of firings (events) of each one.
% GC --> numer of events
% GR --> cell ID

% Important, some cells show no activity in the extinction session, only in retrieval. And vice versa. Therefore, some IDs may not appear.

[GC_R,GR_R] = groupcounts(data.activity{1,2}(:,2));

% In Each cell row
% -->  row: cells
% -->  columns: time

data.activity{2,2} = NaN(max(GR_R),max(GC_R)); % --> time in samples
data.activity{3,2} = NaN(max(GR_R),max(GC_R)); % --> time in seconds
data.activity{4,2} = NaN(max(GR_R),max(GC_R)); % --> event amplitude

data.activity{5,2} = zeros(size(data.activity{2, 2},1),round(data.activity{1, 2}(end,1).*10)); % --> entire session binary format
data.activity{6,2} = zeros(size(data.activity{2, 2},1),round(data.activity{1, 2}(end,1).*10)); % --> entire session binary format

for ii = 1:max(data.activity{1,2}(:,2))

    % Time in samples
    temp_data.activity = data.activity{1,2}(data.activity{1,2}(:,2) == ii,1)';
    data.activity{2,2}(ii,1:length(temp_data.activity)) = round(temp_data.activity * parameters.sampling_rate) ;
    data.activity{2,2}(ii,data.activity{2,2}(ii,:) == 0) = 1; % As the data.activity was transformed from seconds to samples, times of 0 were converted to sample 1.

    % Time in seconds
    data.activity{3,2}(ii,1:length(temp_data.activity)) = temp_data.activity;

    % Binary events over time
    temp_data.activity((temp_data.activity * parameters.sampling_rate) == 0) = 1;     % As the data.activity was transformed from seconds to samples, times of 0 were converted to sample 1.
    data.activity{5,2}(ii,round(temp_data.activity * parameters.sampling_rate)) = 1;  % time sample to binary

    % Amplitude of the event
    temp_data.activity_amplitude = data.activity{1,2}(data.activity{1,2}(:,2) == ii,3)';
    data.activity{4,2}(ii,1:length(temp_data.activity_amplitude)) = temp_data.activity_amplitude;

    %  Amplitude of the event over time
    data.activity{6,1}(ii,round(temp_data.activity * parameters.sampling_rate)) = temp_data.activity_amplitude;  % time sample to binary

end

% Time vector in seconds full session
data.activity{7,2} = linspace(0,data.activity{1, 2}(end,1),round(data.activity{1, 2}(end,1).*parameters.sampling_rate));



clear('GR_E','GR_R','GC_E','GC_R','ii','temp_data')


%% test pca

[idx c sumd] = kmeans(data.activity{1,1}(3:end,[1 2]),10);
[ ~,score] = pca(data.activity{1,1}(3:end,[2 3]));

pc1 = score(:,1)
pc2 = score(:,2)

figure
gscatter(pc1,pc2,idx)

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








