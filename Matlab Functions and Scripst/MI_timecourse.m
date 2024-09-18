
%% Phase-amplitude Cross-frequency coupling measure
% MI over time

% - Performs analysis with raw and surrogate values

% The code relies on the following functions:
% --> ModIndex.m    - by Adriano Tort, Instituto do Cerebro - Universidade Federal do Rio Grande do Norte
% --> shuffle_esc.m - by Rafal Bogacz, Angela Onslow, May 2010

% See Tort et al, 2010 -> 10.1152/jn.00106.2010 
%     Tort et al, 2008 -> 10.1073/pnas.0810524105

% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  01/2024
% Last update: 01/2024


%% To do
% PLOTS and need to be fixed from line 446

%%
% First extract the data with: Extracting_LFPs_and_events.m
% ... then organize the data with the script: Pre_processing.m

%% Run each session sequentially

%% Hilbert Transform

% MI.data

% 1st column - > phase
% 2nd column - > amplitude
% 3nd column - > time vector

% in each cell -> frequencies cutoff
% - rows          - > channels
% - columns       - > time
% - 3nd dimension - > behavioral events

% Choose frequency bands index according to the Pre_processing.m
% MI.params.freq_idx:
% 2 - longcutoff / 3 - theta / 4 - low heta (3-6)
% 5  - hightheta (6 - 9) / 6 - lowgamma (30 50)  / 7 - highgamma (62 100)


MI.params.freq_idx = [3 7]; %frequency 

% Loop over channels and make hilbert transform
% Initializing with NaN
MI.data = cell(3,size(MI.params.freq_idx,2));

MI.data(1,:) = {nan(size(data.lfp{6,1}))}; % baseline
MI.data(2,:) = {nan(size(data.lfp{7,1}))}; % CS-Trials
MI.data(3,:) = {nan(size(data.lfp{8,1}))}; % ITI

% baseline
for ii = 1:length(MI.params.freq_idx)
    for jj = 1:size(data.lfp{6, MI.params.freq_idx(ii)},3)
        for ll = 1:size(data.lfp{6, MI.params.freq_idx(ii)},1)
        
            temp = data.lfp{6, MI.params.freq_idx(ii)}(ll,:,jj);
            temp(isnan(temp(:,:)))=[]; % it was necessary to remove the nan for the function to work

            if ii == 1
                MI.data{1, ii}(ll,1:length(temp),jj) = angle(hilbert(temp));
            else
                MI.data{1, ii}(ll,1:length(temp),jj) = abs(hilbert(temp));
            end
    
        end 
    end 
end


% CS-trials
for ii = 1:length(MI.params.freq_idx)
    for jj = 1:size(data.lfp{7, MI.params.freq_idx(ii)},3)
        for ll = 1:size(data.lfp{7, MI.params.freq_idx(ii)},1)
        
            temp = data.lfp{7, MI.params.freq_idx(ii)}(ll,:,jj);
            temp(isnan(temp(:,:)))=[]; % it was necessary to remove the nan for the function to work

            if ii == 1
                MI.data{2, ii}(ll,1:length(temp),jj) = angle(hilbert(temp));
            else
                MI.data{2, ii}(ll,1:length(temp),jj) = abs(hilbert(temp));
            end
    
        end 
    end 
end


% ITI
for ii = 1:length(MI.params.freq_idx)
    for jj = 1:size(data.lfp{8, MI.params.freq_idx(ii)},3)
        for ll = 1:size(data.lfp{8, MI.params.freq_idx(ii)},1)
        
            temp = data.lfp{8, MI.params.freq_idx(ii)}(ll,:,jj);
            temp(isnan(temp(:,:)))=[]; % it was necessary to remove the nan for the function to work

            if ii == 1
                MI.data{3, ii}(ll,1:length(temp),jj) = angle(hilbert(temp));
            else
                MI.data{3, ii}(ll,1:length(temp),jj) = abs(hilbert(temp));
            end
    
        end 
    end 
end


% Time vector
MI.data{1,3} = linspace(0,size(MI.data{1,1},2)/parameters.decimated_srate,size(MI.data{1, 1},2));
MI.data{2,3} = linspace(0,size(MI.data{2,1},2)/parameters.decimated_srate,size(MI.data{2, 1},2));
MI.data{3,3} = linspace(0,size(MI.data{3,1},2)/parameters.decimated_srate,size(MI.data{3, 1},2));


clear('ii','jj','ll','temp')

%% - CHANNELS MAP - 
% * just in case to check

% For Tucge's data: all channels are in the RE
% For Totty's data: mPFC PL - 01 -> 07
%                   mPFC IL - 08 -> 14
%                   HPC     - 15 -> 29

%% Phase-amplitude cross-frequency coupling measure
 % Continuous time course with overlap
 
% MI_value_timewin 
% - cell          - > first line baseline / second line: CS-Trials / third line: ITI
%                     columns: behavioral events
% - each cell     - > lines: channels
%                     columns: time 

% MeanAmp 
% - cell          - > first line baseline / second line: CS-Trials / third line: ITI
%                     columns: behavioral events
% - each cell     - > lines: channels
%                     columns: amplitude mean values
% - 3nd dimension - > time 


% Time window
MI.params.time_window     =  2; % sec.
MI.params.time_window_idx = round(MI.params.time_window * parameters.decimated_srate);

% Overlap
MI.params.timeoverlap    = .9; % percentage
overlap = round((MI.params.time_window_idx)-(MI.params.timeoverlap * MI.params.time_window_idx));

% Baseline
% Define idx without NaN`s length

MI.params.time2save_idx_baseline = cell(size(MI.data{1,1},3),1);
for ii = 1:size(MI.data{1,1},3)    
    temp_baseline = MI.data{1,2}(:,:,ii);    
    nonantemp_baseline = reshape(temp_baseline(~isnan(temp_baseline)),size(MI.data{1,1},1),[]);
    MI.params.time2save_idx_baseline{ii,1} = (1:overlap:length(nonantemp_baseline) - MI.params.time_window_idx);
end

% CS-Trials
% Define idx without NaN`s length

MI.params.time2save_idx_CS_Trials = cell(size(MI.data{2,1},3),1);
for ii = 1:size(MI.data{2,1},3)    
    temp_CS_Trials = MI.data{2,1}(:,:,ii);    
    nonantemp_CS_Trials = reshape(temp_CS_Trials(~isnan(temp_CS_Trials)),size(MI.data{2,1},1),[]);
    MI.params.time2save_idx_CS_Trials{ii,1} = (1:overlap:length(nonantemp_CS_Trials) - MI.params.time_window_idx);
end

% ITI
% Define idx without NaN`s length

MI.params.time2save_idx_ITI = cell(size(MI.data{3,1},3),1);
for ii = 1:size(MI.data{3,1},3)    
    temp_ITI = MI.data{3,1}(:,:,ii);    
    nonantemp_ITI = reshape(temp_ITI(~isnan(temp_ITI)),size(MI.data{3,1},1),[]);
    MI.params.time2save_idx_ITI{ii,1} = (1:overlap:length(nonantemp_ITI) - MI.params.time_window_idx);
end


% Define number number of phase bins 
MI.params.nbin = 20; 

% variable not centered in the phase bin (rad)
MI.params.position = zeros(1,MI.params.nbin); 

MI.params.winsize = 2*pi/MI.params.nbin;

for jj = 1:MI.params.nbin
    MI.params.position(jj) = -pi+(jj-1)*MI.params.winsize;
end



% Modulation Index for each behavior event

% Baseline
for ii = 1:size(MI.data{1,1},1)
    for jj = 1:size(MI.data{1,1},3)
        for ll = 1:length(MI.params.time2save_idx_baseline{jj})
        
        [MI.MI_value_timewin{1,jj}(ii,ll),MI.MeanAmp_timewin{1,jj}(ii,:,ll)] = ModIndex(MI.data{1,1}(ii,MI.params.time2save_idx_baseline{jj}(ll):(MI.params.time2save_idx_baseline{jj}(ll) + MI.params.time_window_idx -1),jj), MI.data{1,2}(ii,MI.params.time2save_idx_baseline{jj}(ll):(MI.params.time2save_idx_baseline{jj}(ll) + MI.params.time_window_idx -1),jj),MI.params.position);

        end
    end
end 

% CS-Trials
for ii = 1:size(MI.data{2,1},1)
    for jj = 1:size(MI.data{2,1},3)
        for ll = 1:length(MI.params.time2save_idx_CS_Trials{jj})
        
        [MI.MI_value_timewin{2,jj}(ii,ll),MI.MeanAmp_timewin{2,jj}(ii,:,ll)] = ModIndex(MI.data{2,1}(ii,MI.params.time2save_idx_CS_Trials{jj}(ll):(MI.params.time2save_idx_CS_Trials{jj}(ll) + MI.params.time_window_idx -1),jj), MI.data{2,2}(ii,MI.params.time2save_idx_CS_Trials{jj}(ll):(MI.params.time2save_idx_CS_Trials{jj}(ll) + MI.params.time_window_idx -1),jj),MI.params.position);

        end
    end
end 

% ITI
for ii = 1:size(MI.data{3,1},1)
    for jj = 1:size(MI.data{3,1},3)
        for ll = 1:length(MI.params.time2save_idx_ITI{jj})
        
        [MI.MI_value_timewin{3,jj}(ii,ll),MI.MeanAmp_timewin{3,jj}(ii,:,ll)] = ModIndex(MI.data{3,1}(ii,MI.params.time2save_idx_ITI{jj}(ll):(MI.params.time2save_idx_ITI{jj}(ll) + MI.params.time_window_idx -1),jj), MI.data{3,2}(ii,MI.params.time2save_idx_ITI{jj}(ll):(MI.params.time2save_idx_ITI{jj}(ll) + MI.params.time_window_idx -1),jj),MI.params.position);

        end
    end
end 

% Time vector
MI.MI_timev_timewin{1,1} = linspace(0,size(MI.data{1,1},2)/parameters.decimated_srate,size(MI.MI_value_timewin{1,1},2));
MI.MI_timev_timewin{2,1} = linspace(0,size(MI.data{2,1},2)/parameters.decimated_srate,size(MI.MI_value_timewin{2,1},2));
MI.MI_timev_timewin{3,1} = linspace(0,size(MI.data{3,1},2)/parameters.decimated_srate,size(MI.MI_value_timewin{3,1},2));

clear ('overlap','nonan','temp','MI.params.time2save_idx_1','MI.params.time2save_idx_2','ii','jj','ll','temp1','nonan_len')

%% Plot MI values over time for each behaviour event

% Choose channel to plot
ch = 14;

% Set Figure
figure
set(gcf,'color','white')

subplot(2,size(MI.MI_value_timewin,2)+1,1)
sgtitle({'Phase-Amplitude coupling';['(Sliding window = ' num2str(MI.params.time_window) 's' ' - ' 'overlap = ' num2str(MI.params.timeoverlap*100) '%)']}) 

plot(MI.MI_timev_timewin{1,1}, MI.MI_value_timewin{1,1}(ch,:), 'Color','[0.6350, 0.0780, 0.1840]','linew',1)
%axis square
set(gca,'fontsize',11)
xlabel('Time (s)')
ylabel('Modulation Index')
%xlim([])
ylim([0 0.015])

title({'Baseline';[]})

box off


for ii = 1:size(MI.MI_value_timewin,2)
    subplot(2,size(MI.MI_value_timewin,2)+1,ii+1)
    sgtitle({[];'Phase-Amplitude coupling';['(Sliding window = ' num2str(MI.params.time_window) 's' ' - ' 'overlap = ' num2str(MI.params.timeoverlap*100) '%)']})

    plot(MI.MI_timev_timewin{2,1}, MI.MI_value_timewin{2,ii}(ch,:), 'Color','[0.6350, 0.0780, 0.1840]','linew',1)
    %axis square
    set(gca,'fontsize',11)
    xlabel('Time (s)')
    ylabel('Modulation Index')
    %xlim([])
    ylim([0 0.015])

    title ({['CS: ' num2str(ii)];[]})

    box off
end

for ii = 1:size(MI.MI_value_timewin,2)
    subplot(2,size(MI.MI_value_timewin,2)+1,ii+7)
    sgtitle({[];'Phase-Amplitude coupling';['(Sliding window = ' num2str(MI.params.time_window) 's' ' - ' 'overlap = ' num2str(MI.params.timeoverlap*100) '%)']})

    plot(MI.MI_timev_timewin{3,1}, MI.MI_value_timewin{3,ii}(ch,:), 'Color','[0.6350, 0.0780, 0.1840]','linew',1)
    %axis square
    set(gca,'fontsize',11)
    xlabel('Time (s)')
    ylabel('Modulation Index')
    %xlim([])
    ylim([0 0.015])

    title ({['ITI: ' num2str(ii)];[]})

    box off
end



clear ('ii','ch')

%% Surrogate phase vectors to compare

% MI_shuffle_values_timewin
% - cell          - > first line: baseline / second line: CS-Trials / third line: ITI
%                     column: behavioral events
% - each cell     - > lines: channels
%                     columns: time 
%                     3nd dimension - > rearrangements - each timestamp idx

numshf  = 200; % number of shuffled segments
nsurrog = 200; % number of rearrangements

% Loop over events 

baseline_shuffle = [];
CS_shuffle = [];


% Baseline
for ss = 1:nsurrog
    
    for ii = 1:size(MI.data{1,1},1)
        for jj = 1:size(MI.data{1,1},3)
            for ll = 1:length(MI.params.time2save_idx_baseline{1,1})
        
            baseline_shuffle{1,jj}(ii,:,ll) = shuffle_esc(MI.data{1,1}(ii,MI.params.time2save_idx_baseline{jj}(ll):(MI.params.time2save_idx_baseline{jj}(ll) + MI.params.time_window_idx -1),jj),parameters.decimated_srate,numshf);
            
            end
        end
    end
    
    for ii = 1:size(MI.data{1,1},1)
        for jj = 1:size(MI.data{1,1},3)
            for ll = 1:length(MI.params.time2save_idx_baseline{1,1})
                
            [MI.MI_shuffle_values_timewin{1,jj}(ii,ll,ss)] = ModIndex(baseline_shuffle{1,jj}(ii,:,ll), MI.data{1,2}(ii,MI.params.time2save_idx_baseline{jj}(ll):(MI.params.time2save_idx_baseline{jj}(ll) + MI.params.time_window_idx -1),jj),MI.params.position);

            end
        end
    end
    
        baseline_shuffle = [];
        
end 


% CS_Trials
for ss = 1:nsurrog
    
    for ii = 1:size(MI.data{2,1},1)
        for jj = 1:size(MI.data{2,1},3)
            for ll = 1:length(MI.params.time2save_idx_CS_Trials{1,1})
        
            CS_Trials_shuffle{1,jj}(ii,:,ll) = shuffle_esc(MI.data{2,1}(ii,MI.params.time2save_idx_CS_Trials{jj}(ll):(MI.params.time2save_idx_CS_Trials{jj}(ll) + MI.params.time_window_idx -1),jj),parameters.decimated_srate,numshf);
            
            end
        end
    end
    
    for ii = 1:size(MI.data{2,1},1)
        for jj = 1:size(MI.data{2,1},3)
            for ll = 1:length(MI.params.time2save_idx_CS_Trials{1,1})
                
            [MI.MI_shuffle_values_timewin{2,jj}(ii,ll,ss)] = ModIndex(CS_Trials_shuffle{1,jj}(ii,:,ll), MI.data{2,2}(ii,MI.params.time2save_idx_CS_Trials{jj}(ll):(MI.params.time2save_idx_CS_Trials{jj}(ll) + MI.params.time_window_idx -1),jj),MI.params.position);

            end
        end
    end
    
        CS_Trials_shuffle = [];
        
end 


% ITI
for ss = 1:nsurrog
    
    for ii = 1:size(MI.data{3,1},1)
        for jj = 1:size(MI.data{3,1},3)
            for ll = 1:length(MI.params.time2save_idx_ITI{1,1})
        
            ITI_shuffle{1,jj}(ii,:,ll) = shuffle_esc(MI.data{3,1}(ii,MI.params.time2save_idx_ITI{jj}(ll):(MI.params.time2save_idx_ITI{jj}(ll) + MI.params.time_window_idx -1),jj),parameters.decimated_srate,numshf);
            
            end
        end
    end
    
    for ii = 1:size(MI.data{3,1},1)
        for jj = 1:size(MI.data{3,1},3)
            for ll = 1:length(MI.params.time2save_idx_ITI{1,1})
                
            [MI.MI_shuffle_values_timewin{3,jj}(ii,ll,ss)] = ModIndex(ITI_shuffle{1,jj}(ii,:,ll), MI.data{3,2}(ii,MI.params.time2save_idx_ITI{jj}(ll):(MI.params.time2save_idx_ITI{jj}(ll) + MI.params.time_window_idx -1),jj),MI.params.position);

            end
        end
    end
    
        ITI_shuffle = [];
        
end 

% z surrogated values
% baseline    
MI.stats.z_MI_shuffle_values_timewin{1,1} = zscore(MI.MI_shuffle_values_timewin{1,1}); % baseline

for ii = 1:size(MI.MI_shuffle_values_timewin,2)
    MI.stats.z_MI_shuffle_values_timewin{2,ii} = zscore(MI.MI_shuffle_values_timewin{2,ii}); % CS-Trials
    MI.stats.z_MI_shuffle_values_timewin{3,ii} = zscore(MI.MI_shuffle_values_timewin{3,ii}); % ITI
end


% z real(observed) values
% baseline
MI.stats.z_MI_value_timewin{1,1} = (MI.MI_value_timewin{1,1} - mean(MI.MI_shuffle_values_timewin{1,1},3))./std(MI.MI_shuffle_values_timewin{1,1},[],3); % baseline

for ii = 1:size(MI.MI_shuffle_values_timewin,2)
    MI.stats.z_MI_value_timewin{2,ii} = (MI.MI_value_timewin{2,ii} - mean(MI.MI_shuffle_values_timewin{2,ii},3))./std(MI.MI_shuffle_values_timewin{2,ii},[],3); % CS-Trials
    MI.stats.z_MI_value_timewin{3,ii} = (MI.MI_value_timewin{3,ii} - mean(MI.MI_shuffle_values_timewin{3,ii},3))./std(MI.MI_shuffle_values_timewin{3,ii},[],3); % ITI
end


% p values
% baseline
MI.stats.p_MI_value_timewin{1,1} = 2*(1-normcdf(MI.stats.z_MI_value_timewin{1,1}));

for ii = 1:size(MI.MI_shuffle_values_timewin,2)     
    MI.stats.p_MI_value_timewin{2,ii} = 2*(1-normcdf(MI.stats.z_MI_value_timewin{2,ii}));
    MI.stats.p_MI_value_timewin{3,ii} = 2*(1-normcdf(MI.stats.z_MI_value_timewin{3,ii}));

%    zval = norminv(1-(.05/length(sessions)));                                                      % z-value threshold at p=0.05, correcting for multiple comparisons
%    MI.stats.p_MI_value_timewin{1,1} = 0.5 * erfc(-MI.stats.z_MI_value_timewin{1,1} ./ sqrt(2));   % p value. Similar to Matlab function: normcdf(-z) two-tailed
%    MI.stats.p_MI_value_timewin{2,1} = 0.5 * erfc(-MI.stats.z_MI_value_timewin{2,ii} ./ sqrt(2));  % p value. Similar to Matlab function: normcdf(-z) two-tailed
%    MI.stats.p_MI_value_timewin{3,1} = 0.5 * erfc(-MI.stats.z_MI_value_timewin{3,ii} ./ sqrt(2));  % p value. Similar to Matlab function: normcdf(-z) two-tailed

end

clear('ll', 'ss','jj','ii','nsurrog','numshf','time2samples','baseline_shuffle','CS_shuffle')

%% Plot MI values over time for each behaviour event

% Choose channel to plot
ch = 6;

% identify p values idx

p = 0.05;

p_idx{1,1} = MI.stats.p_MI_value_timewin{1, 1} < p;

for ii = 1:size(MI.stats.p_MI_value_timewin,2)   
    p_idx{2,ii} = MI.stats.p_MI_value_timewin{2, ii} < p;
    p_idx{3,ii} = MI.stats.p_MI_value_timewin{3, ii} < p;    
end

% Set Figure
figure
set(gcf,'color','white')

suptitle({[];'Phase-Amplitude coupling';['(Sliding window = ' num2str(MI.params.time_window) 's' ' - ' 'overlap = ' num2str(MI.params.timeoverlap*100) '%)']}) 

% Plot RAW MI values
for ii= 1:size(MI.MI_value_timewin,1)
    
    subplot(2,size(MI.MI_value_timewin,1),ii)
    
    t1   = [MI.MI_value_timewin_time_vector{ii,1} MI.MI_value_timewin_time_vector{ii,2}];
    plt1 = [MI.MI_value_timewin{ii,1}(ch,:) MI.MI_value_timewin{ii,2}(ch,:)];
    
    % Concatenate period to plot
    plot(t1,plt1, 'Color','k','linew',2)
    hold all
    plot(t1([p_idx{ii,1}(ch,:) p_idx{ii,2}(ch,:)]), plt1([p_idx{ii,1}(ch,:) p_idx{ii,2}(ch,:)]), 'Ro','linew',2)    
    plot([0 0],[-5 5],'k--')
    ylim([0 0.03])
    
    title(['Bout ',num2str(ii)])
    xlabel('Time (s)')
    ylabel('Modulation Index')
    axis square
    box off

end


% Plot MI z values
for ii= 1:size(MI.MI_value_timewin,1)
    
    subplot(2,size(MI.MI_value_timewin,1),ii+size(MI.MI_value_timewin,1))
    
    t2   = [MI.MI_value_timewin_time_vector{ii,1} MI.MI_value_timewin_time_vector{ii,2}];
    plt2 = [MI.stats.z_MI_value_timewin{ii,1}(ch,:) MI.stats.z_MI_value_timewin{ii,2}(ch,:)];
    
    % Concatenate period to plot
    plot(t2,plt2, 'Color','k','linew',2)
    hold all
    plot(t2([p_idx{ii,1}(ch,:) p_idx{ii,2}(ch,:)]), plt2([p_idx{ii,1}(ch,:) p_idx{ii,2}(ch,:)]), 'Ro','linew',2)    
    plot([t2(1) t2(end)],[1.96 1.96],'k--')
    plot([t2(1) t2(end)],[-1.96 -1.96],'k--')
    plot([0 0],[-5 5],'k--')
    ylim([-6 6])
    
    title(['Bout ',num2str(ii)])
    xlabel('Time (s)')
    ylabel('Z values')
    axis square
    box off

end

    legend('Theta (4–10 Hz) & Gamma (80–140 Hz)','location','southoutside')
    legend boxoff

clear ('ii','ch','g','idx','p','p_idx','t1','plt1','t2','plt2')


%% Plot one trial and choose desire window amplitude distribution

% Choose trial
trial = 1;

% Choose channel
ch = 3;

% Concatenate data and time 
time = [MI.MI_value_timewin_time_vector{trial,1} MI.MI_value_timewin_time_vector{trial,2}];

% z values
data_z = [MI.stats.z_MI_value_timewin{trial,1}(ch,:) MI.stats.z_MI_value_timewin{trial,2}(ch,:)];

% MI raw values
data_mi = [MI.MI_value_timewin{trial,1}(ch,:) MI.MI_value_timewin{trial,2}(ch,:)];


% Plot to choose best time windows. 
% Use data_cursor and export cursor_info values to the workspace 
figure
set(gcf,'color','white')

subplot(2,1,1)
plot(time,data_mi, 'Color','k','linew',2)
hold
plot([0 0],[0 max(data_mi)],'k--')

subplot(2,1,2)
plot(time,data_z, 'Color','k','linew',2)
hold
plot([time(1) time(end)],[1.96 1.96],'k--')
plot([time(1) time(end)],[-1.96 -1.96],'k--')
plot([0 0],[min(data_z) max(data_z)],'k--')

%Define IDX
idx = fliplr([cursor_info.DataIndex]);

% Concatenate amplitude histogram
% - each column is a differente time point
hstg_1 = squeeze(cat(3,MI.MeanAmp_timewin{trial,1}(ch,:,:), MI.MeanAmp_timewin{trial,2}(ch,:,:)));

% Concatenate MI shuffle z values distribuition 
% - each line is a differente time point
hstg_2 = squeeze(cat(2,MI.stats.z_MI_shuffle_values_timewin{trial,1}(ch,:,:), MI.stats.z_MI_shuffle_values_timewin{trial,2}(ch,:,:)));

% Concatenate MI real z values 
% - each column is a differente time point
MI_z = squeeze(cat(2,MI.stats.z_MI_value_timewin{trial,1}(ch,:), MI.stats.z_MI_value_timewin{trial,2}(ch,:)));

% Define x values for amplitude histogram
xvalue1 = rad2deg(MI.params.position) + 180;
xvalue2 = xvalue1 + 360;


figure (2)
set(gcf,'color','white')

suptitle({[];'Phase-Amplitude coupling';['(Sliding window = ' num2str(MI.params.time_window) 's' ' - ' 'overlap = ' num2str(MI.params.timeoverlap*100) '%)']}) 

subplot(4,length(idx),[1 length(idx)])
plot(time,data_mi, 'Color','k','linew',2)
hold
plot(time(idx),data_mi(idx),'Ro','linew',4,'MarkerFaceColor','R')
plot([0 0],[min(data_mi) max(data_mi)],'k--')

title(['\fontsize{12}Bout' num2str(trial)])
xlabel('\fontsize{12}Time (s)')
ylabel('\fontsize{12}MI')
box off

subplot(4,length(idx),[length(idx)+1 length(idx)*2])
plot(time,data_z, 'Color','k','linew',2)
hold
plot([time(1) time(end)],[1.96 1.96],'k--')
plot([time(1) time(end)],[-1.96 -1.96],'k--')

plot([0 0],[min(data_z) max(data_z)],'k--')
plot(time(idx),data_z(idx),'Ro','linew',4,'MarkerFaceColor','R')
xlabel('\fontsize{12}Time (s)')
ylabel('\fontsize{12}z-value')
ylim([-5 5])
box off

for ii = 1:length(idx)
    subplot(4,length(idx),2*length(idx)+ii)
    b3 = bar([xvalue1 xvalue2],[hstg_1(:,idx(ii)); hstg_1(:,idx(ii))]);
    b3.FaceColor = 'w';
    ylim([0 0.1])
    xlabel('\fontsize{12}phase (degree)')
    ylabel('\fontsize{12}norm amplitude')
end

for ii = 1:length(idx)
    subplot(4,length(idx),3*length(idx)+ii)
    h = histogram(hstg_2(ii,:),MI.params.nbin);
    h.FaceColor = 'w';
    %xlim([-3 +3])
    hold on
    plot([MI_z(1,idx(ii)) MI_z(1,idx(ii))],[0 30],'Color', '[0.6350, 0.0780, 0.1840]','linewidth',4);
    xlabel('\fontsize{12}z values')
    ylabel('\fontsize{14}Frequency')
    legend({'\fontsize{12}Random';'\fontsize{12}Observed'},'box','off');
end

clear('ch','data_mi','data_z','time','trial','cursor_info','idx','hstg_1','hstg_2','MI_z','xvalue1','xvalue2','b3','h','ii')


%% Plot all channels over all trials. 

% MI raw values
figure
suptitle({[];'Phase-Amplitude coupling';['(Sliding window = ' num2str(MI.params.time_window) 's' ' - ' 'overlap = ' num2str(MI.params.timeoverlap*100) '%)']}) 


for ii= 1:size(MI.MI_value_timewin,1)    
    subplot(1,size(MI.MI_value_timewin,1),ii)
    mi_values = [MI.MI_value_timewin{ii, 1} MI.MI_value_timewin{ii, 2}];
    t = [MI.MI_value_timewin_time_vector{ii,1} MI.MI_value_timewin_time_vector{ii,2}];
    
    contourf(t,(1:1:size(mi_values,1)),mi_values,80,'linecolor','none');
    c = colorbar;
    ylabel(c,'MI','FontSize',12,'Rotation',270);
    c.Label.Position = [4.88,0.015,0];
    caxis([0 0.03])
    xlabel('\fontsize{12}time (s)')
    ylabel('\fontsize{12}channels')

end


% p values
figure

for ii= 1:size(MI.MI_value_timewin,1)
    
    subplot(1,size(MI.MI_value_timewin,1),ii)
    p_values = [MI.stats.p_MI_value_timewin{ii, 1} MI.stats.p_MI_value_timewin{ii, 2}];
    t = [MI.MI_value_timewin_time_vector{ii,1} MI.MI_value_timewin_time_vector{ii,2}];

    contourf(t,(1:1:size(p_values,1)),p_values,80,'linecolor','none');
    colormap(flipud(jet))
    c = colorbar;
    ylabel(c,'p values','FontSize',12,'Rotation',270);
    c.Label.Position = [4.88,0.049,0];
    caxis([0 0.1])
    xlabel('\fontsize{12}time (s)')
    ylabel('\fontsize{12}channels')
    
end

clear('ii','c','mi_values','p_values','t')

%% Video MI

% Choose channel
MI.video.ch = 3;

% choose trial
MI.video.trl = 1;

% concatenate period

% MI raw values
MI.video.data_2plot = [MI.MI_value_timewin{MI.video.trl,1}(MI.video.ch,:) MI.MI_value_timewin{MI.video.trl,2}(MI.video.ch,:)];
% MI z values
%MI.video.data_2plot = [MI.stats.z_MI_value_timewin{MI.video.trl,1}(MI.video.ch,:) MI.stats.z_MI_value_timewin{MI.video.trl,2}(MI.video.ch,:)];
% Time vector
MI.video.time_2plot = [MI.MI_value_timewin_time_vector{MI.video.trl,1} MI.MI_value_timewin_time_vector{MI.video.trl,2}];

% set video time idx
MI.video.vid_min = -5; % check pre trial time defined before
MI.video.vid_min_idx = dsearchn(MI.video.time_2plot',MI.video.vid_min);

MI.video.vid_max     = 20; % check trial time defined before
MI.video.vid_max_idx = dsearchn(MI.video.time_2plot',MI.video.vid_max);

% Define Plot (parameters for high resolution)
figure('units', 'pixels', 'position', [0 0 1920 1080]); clf
set(gcf,'color','white')
axis([-5,20,0,0.025])

set(gca,'XColor','w','Fontsize',14)

%xlabel('Time (s)')
ylabel('Modulation Index')
axis square
hold all
box off

plot([18 20],[0.0005 0.0005],'k','linew',2)
plot([0 0],[0 0.025],'k--','linew',.5)

h = animatedline;

h.Color = [0.6350, 0.0780, 0.1840];
h.LineWidth = 3;

% setup movie
mov = VideoWriter('MI','Uncompressed AVI');

open(mov)

tic
for ii = 1:length(MI.video.vid_min_idx:MI.video.vid_max_idx)
    
    addpoints(h,MI.video.time_2plot(ii),MI.video.data_2plot(ii));
    drawnow
    pause (0.220)
    
    F = getframe(gcf);
    writeVideo(mov,F);
end
toc


close(mov)

clear('F','h','ii','mov')

%% Video amplitude mean values

% Channel and trial define before in Video MI over time

% Concatenate period
MI.video.MeanAmp_2plot = squeeze(cat(3,MI.MeanAmp_timewin{MI.video.trl,1}(MI.video.ch,:,:), MI.MeanAmp_timewin{MI.video.trl,2}(MI.video.ch,:,:)));

% Define Plot (parameters for high resolution)
figure('units', 'pixels', 'position', [0 0 1920 1080]); clf
set(gcf,'color','w');

% Set x axis
xvalue1 = (rad2deg(MI.params.position) + 180)';
xvalue2 = (xvalue1 + 360);

% Define pre values as zeros
b = bar(cat(1, xvalue1, xvalue2),zeros(length(xvalue1) * 2,1));
b.FaceColor = [0.8, 0.8, 0.8];
axis([0,xvalue2(end),0,0.15])

ylabel('\fontsize{14} Amplitude')
xlabel('\fontsize{14} Phase')
set(gca,'Fontsize',14)

axis square
box off

mov = VideoWriter('amp','Uncompressed AVI');
open(mov)

tic
for  ii = 1:size(MI.video.MeanAmp_2plot,2)
     
     % catch values
     amp = cat(1,MI.video.MeanAmp_2plot(:,ii), MI.video.MeanAmp_2plot(:,ii));
     
     % set values in each iteneration
     set(b,'YData',amp)
     
     pause(0.110)
     
     F = getframe(gcf);
     writeVideo(mov,F);
      
end
toc

close(mov)

clear('amp','xvalue1','xvalue2','b','F','h','ii','mov')


%% last update 20/09/2020 - 22:16
%  listening: Grouper - Headache
