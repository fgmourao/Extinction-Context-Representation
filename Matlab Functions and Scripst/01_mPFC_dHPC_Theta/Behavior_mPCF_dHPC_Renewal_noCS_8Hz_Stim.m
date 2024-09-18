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
% Last update: 05/2024

%% Trial epochs - baseline and CS sound period
% Original sample rate

% Freezing settings
parameters.thr_1 = 5; % lower threshold in percentage.
parameters.thr_2 = 1; % higher threshold in sec to consider freezing.

% Full session Freezing

% idx for freezing epochs
[data.behavior{2,1}(1,:), data.behavior{2,1}(2,:), ~] = ZeroOnesCount(data.behavior{1,1} <= parameters.thr_1); % first row = start events / second = time in samples
idx_to_remove = data.behavior{2, 1}(2,:) < parameters.thr_2 * parameters.original_srate; % Exclude values based on timing threshold
data.behavior{2, 1}(:,idx_to_remove) = [];                                      
data.behavior{2, 1}(3,:) = data.behavior{2, 1}(2,:)./parameters.original_srate; % thirt row = time in second


% idx for non-freezing epochs
% First lets find all freezing idx
f_idx(:,1) = data.behavior{2, 1}(1,:);
f_idx(:,2) = data.behavior{2, 1}(1,:) + data.behavior{2, 1}(2,:);

%Then, lets find the intervals who represents non-freezing
data.behavior{2,2}(:,1) = [1 f_idx(1,1)-1]; % from beggining until the first freezing event

if size(f_idx,1)>1
    data.behavior{2,2}(1,2:length(f_idx)) = f_idx(1:end-1,2)+1;                      % first row = start events
    data.behavior{2,2}(2,2:length(f_idx)) = (f_idx(2:end,1)-1) - f_idx(1:end-1,2)-1; % second = time in samples
end

data.behavior{2,2}(3,:) = data.behavior{2, 2}(2,:)./parameters.original_srate;   % thirt row = time in second

idx_to_remove_1 = data.behavior{2, 2}(2,:) < parameters.thr_2 * parameters.original_srate; % Exclude values based on timing threshold
data.behavior{2, 2}(:,idx_to_remove_1) = []; 

clear('f_idx','idx_to_remove','idx_to_remove_1')

%% Trial epochs

% baseline
data.behavior{3,1}  = data.behavior{1,1}(1,1:data.events{2, 1}(1,1)-1);

% CS sound and ITI period
for ii = 1:size(data.events{2, 1},1)

        epochs{1,ii}  = data.behavior{1,1}(1,data.events{2, 1}(ii,1):data.events{2, 1}(ii,2)); % CS-TONES
        epochs{2,ii}  = data.behavior{1,1}(1,data.events{2, 2}(ii,1):data.events{2, 2}(ii,2)); % ITI-Trials

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

Baseline_time = 180;

% Freezing idx
for ii = 1:size(data.behavior{2,1},2)

    if data.behavior{2,1}(1,ii) < Baseline_time * parameters.original_srate; % Do not use baseline idx
        continue
    end

    data.events_behavior{1,1}(ii,1) = data.behavior{2,1}(1,ii); % Start freezezing
    data.events_behavior{1,1}(ii,2) = data.behavior{2,1}(1,ii) + data.behavior{2,1}(2,ii); % End freezezing

end

% delete zeros rows
data.events_behavior{1, 1}( ~any(data.events_behavior{1, 1},2), : ) = [];  %rows

% non-Freezing idx
data.events_behavior{1,2}(:,2) = data.events_behavior{1,1}(1:end-1,2)+1;
data.events_behavior{1,2}(:,1) = data.events_behavior{1,1}(2:end,1)-1;

% Exclude time epochs < 1s
idx_to_remove = (data.events_behavior{1,2}(:,2) - data.events_behavior{1,2}(:,1)) < 1*parameters.original_srate;
data.events_behavior{1,2}(idx_to_remove,:) = [];

clear('Baseline_time','ii','idx_to_remove')

%% Trial epochs - Binarized data in time epochs
% Binarized data in time epochs just as analyzed in the laboratory.

% -------------------------------------------
% Important Note

% At Plexon, Tugce delete some samples (before and after recording)...to fit? I don`t know why... 
% and exports data with bins of 200 ms (averaging 200 ms time epoch from their original sampling). 
% This represents an output similar to the conditioning box (5Hz sampling rate).
% The data in its original sampling frequency recorded in Plexon versus "Binarized" by averaging, 
% may express slightly differently.

% -------------------------------------------

% Averaging time epochs at each 200ms(5hz)
parameters.time_behav_bins = .2; %sec

data_bins = 1:parameters.time_behav_bins * parameters.original_srate:length(data.behavior{1,1});

for ii = 1:length(data_bins)-1
    data.behavior_bins{1,1}(1,ii) = mean(data.behavior{1,1}(1,data_bins(ii):data_bins(ii+1)-1));
end



% Freezing settings
parameters.thr_1 = 5; % lower threshold in percentage.
parameters.thr_2 = 1; % higher threshold in sec to consider freezing.

% Full session Freezing

[data.behavior_bins{2,1}(1,:), data.behavior_bins{2,1}(2,:), ~] = ZeroOnesCount(data.behavior_bins{1,1} <= parameters.thr_1);

idx_to_remove = data.behavior_bins{2, 1}(2,:) < (1/parameters.time_behav_bins);
data.behavior_bins{2, 1}(:,idx_to_remove) = [];
data.behavior_bins{2, 1}(3,:) = data.behavior_bins{2, 1}(2,:) ./ (1/parameters.time_behav_bins);

clear('ii','idx_to_remove')

%% Trial epochs

% baseline
data.behavior_bins{3,1}  = data.behavior_bins{1,1}(1,1:(data.events{2, 1}(1,1)-1)/(parameters.time_behav_bins * parameters.original_srate));

% CS sound and ITI period
for ii = 1:size(data.events{2, 1},1)

        epochs_bins{1,ii}  = data.behavior_bins{1,1}(1,data.events{2, 1}(ii,1)/(parameters.time_behav_bins * parameters.original_srate) : (data.events{2, 1}(ii,2)-1)/(parameters.time_behav_bins * parameters.original_srate)); % CS-Trials
        epochs_bins{2,ii}  = data.behavior_bins{1,1}(1,data.events{2, 2}(ii,1)/(parameters.time_behav_bins * parameters.original_srate) : (data.events{2, 2}(ii,2)-1)/(parameters.time_behav_bins * parameters.original_srate)); % ITI-Trials

end

% Reshaped CS and ITI trials in correct order and add to data behavior epochs
CS_ITI_Trials_bins = reshape(epochs_bins,1,[]);
data.behavior_bins(3,2:length(CS_ITI_Trials_bins)+1) = CS_ITI_Trials_bins;


% CS-Trials and ITI Freezing
for ii = 1:size(data.behavior_bins,2)
    
    data.behavior_bins{4,ii} = data.behavior_bins{3,ii} <= parameters.thr_1;

    [data.behavior_bins{5,ii}(1,:), data.behavior_bins{5,ii}(2,:), ~] = ZeroOnesCount(data.behavior_bins{4,ii});

    idx_to_remove = data.behavior_bins{5, ii}(2,:) < parameters.thr_2 * (1/parameters.time_behav_bins);
    data.behavior_bins{5, ii}(:,idx_to_remove) = [];
    data.behavior_bins{5, ii}(3,:) = data.behavior_bins{5, ii}(2,:) ./ (1/parameters.time_behav_bins);

    idx_to_remove = [];

    data.behavior_bins{6, ii} = sum(data.behavior_bins{5, ii}(3,:));
    data.behavior_bins{7, ii} = (data.behavior_bins{6, ii}.*100) ./ (length(data.behavior_bins{3,ii}) ./ (1/parameters.time_behav_bins));

end

clear('parameters.time_behav_bins','data_bins','ii','epochs_bins','CS_ITI_Trials_bins','parameters.thr_1','parameters.thr_2','idx_to_remove')

%% Select data to plot

% Movement
data_2_plot_1     = data.behavior_bins{1,1};

% Time vector
behav_bins_time_v = linspace(1,length(data.behavior_bins{1,1})/(1/parameters.time_behav_bins),length(data.behavior_bins{1,1}));

% CS indexes
cs_trial       = round([data.events{2, 1}./(parameters.time_behav_bins * parameters.original_srate)]); 

% Freezing indexes
freezing_start = data.behavior_bins{2, 1}(1,:);
freezing_end   = data.behavior_bins{2, 1}(1,:)+data.behavior_bins{2, 1}(2,:)-1;


% Freezing Percentage
data_2_plot_2  = cell2mat(data.behavior_bins(7, :));


figure
set(gcf,'color','white')
sc = [1,1,1300,980];
set(gcf, 'Position', sc);

subplot(1,2,1)
hold all
sgtitle({'Renewal Recordings - no CS-Tones. Only 8Hz RE opto Stimulation';[]})

%plot(behav_bins_time_v,data_2_plot,'Color',[0.3, 0.3, 0.3, 0.3]) % raw data
plot(behav_bins_time_v, movmean(data_2_plot_1,(1/parameters.time_behav_bins)),'linew', 1, 'Color',[0.6, 0.6, 0.6, 0.6]);
plot(behav_bins_time_v,ones(1,length(behav_bins_time_v)).*5,'k--');
plot([behav_bins_time_v(freezing_start);behav_bins_time_v(freezing_end)], [ones(1,length(freezing_start)).*50;ones(1,length(freezing_end)).*50],'k-','linew', 2,'Color',[.6, 0, 0]);

plot([behav_bins_time_v(cs_trial(:,1));behav_bins_time_v(cs_trial(:,2))], [ones(1,length(cs_trial(:,1))).*95;ones(1,length(cs_trial(:,2))).*95],'k-','linew', 20,'Color',[.5, .7, 1]);

%plot(behav_bins_time_v(100:end-100),ones(1,length(behav_bins_time_v(100:end-100))).*95,'k-','linew', .5,'Color',[1, .4, .4])

xlabel('Time (s)','FontSize',14), ylabel('Movement (%)','FontSize',14)
xlim([behav_bins_time_v(1)-10 behav_bins_time_v(end)])
ylim([0 100])

legend('Movement','Threshold','Freezing','NumColumns',3,'Location','southoutside','FontSize',14)
legend('boxoff')


subplot(1,2,2)
plot(data_2_plot_2,'-go','linew', 2, 'Color',[0 0 .6],'MarkerEdgeColor',[0 0 .6],'MarkerFaceColor',[0 0 .6])    

ylabel('Freezing (%)','FontSize',14)
xticks(1:length(data_2_plot_2))

xticklabels({'Baseline','Stim 1','ITI','Stim 2','ITI','Stim 3','ITI','Stim 4','ITI','Stim 5','ITI','Stim 6','ITI','Stim 7','ITI','Stim 8','ITI','Stim 9','ITI','Stim 10','End'})
xlabel('[]','FontSize',5)
title('All Trials','FontSize',14)
xtickangle(90)
xlim([0.2 length(data_2_plot_2)+0.2])
ylim([-5 105])
box off


% Clear
clear('behav_bins_time_v','cs_trial','data_2_plot_1','data_2_plot_2','freezing_end','freezing_start')

%% Save

newStr = id(1:end-5);

path = '/Users/flavio/Desktop';
%path = files.FilesLoaded{1, 1}.folder;

name = strcat(path,'/',newStr,'_behavior');

saveas(gcf,name,'png')

close all

clear('name','newStr1','path')


%% last update 14/05/2024 - 18:54
%  listening: 

