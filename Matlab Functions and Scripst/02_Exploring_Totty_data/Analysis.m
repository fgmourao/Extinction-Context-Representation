%% Analysis

% Script to analyze all animals
% - The code relies on the following functions : -> Open_Files.m
%                                                -> Extracting_LFPs_and_events_from_all.m
%                                                -> ... other desired functions

% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  12/2023
% Last update: 06/2024

%%

tic
set(groot,'DefaultFigureColormap',jet)

[files] = Open_files();

for ms = 1:length(files.FilesLoaded{1, 1}) % Loop over files path

    % Settings
    id          = files.id(ms).name;
    FilesLoaded = {files.FilesLoaded{1,1}(ms).name};
    Path        = files.FilesLoaded{1,1}(ms).folder;

    %% Open and extract raw files
    %
    %     %Load function
     [data,parameters] = Extracting_LFPs_and_events(id,FilesLoaded,Path);

    %% Scripts pre-processing

    %    check_channels          % Check raw data
    %    Pre_processing_Totty;   % Pre processing
    %    Behavior                % Behavior Analyse
    %    Behavior_exposure       % Behavior Analyse forexposure sessions
    %    CorCov                  % Correlation and Covariance Matrices betwwen channels
    %    Main_Plots              % plot to check raw data

    %% Context Exposure 1

    % Load. Choose data from data.lfp according pre_processing.m and main plots script
    %     name = strcat(Path,'/',id);
    %     fprintf('\n loading data... \n');
    %     load(name,'data','parameters');

    %     % Time window
    %     CSIT{1} = 1:10;      % MT1
    %     CSIT{2} = 1:10;      % MT3
    %     CSIT{3} = 1:10;      % MT4
    %     CSIT{4} = 1:10;      % MT5
    %     CSIT{5} = 1:10;      % MT6
    %     CSIT{6} = 1:10;      % MT7

    %     % Scripts

    %% Extinction Session 1

    % Load. Choose data from data.lfp according pre_processing.m and main plots script
%             name = strcat(Path,'/',id);
%             fprintf('\n loading data... \n');
%
%            load(name,'data','parameters');
%             load(name);


    % Baseline
    % This follows the order of animals according to the main loop over ms values (see Analysis.m). Each session needs to be double-checked.

    % - IMPORTANT -
    % The baseline is being analyzed without noise or freezing periods.
    % All CS and ITI are being saved entirely. Only the plots are being saved without noise trials.


%     % 2) Baseline timeepochs -> first minute. All animals checked. No noise
%             B_clean{1} = dsearchn(data.timev_decimated',[10 170]');   % MT1
%             B_clean{2} = dsearchn(data.timev_decimated',[10 170]');   % MT3
%             B_clean{3} = dsearchn(data.timev_decimated',[10 170]');   % MT4
%             B_clean{4} = dsearchn(data.timev_decimated',[10 170]');   % MT5
%             B_clean{5} = dsearchn(data.timev_decimated',[10 170]');   % MT6
%             B_clean{6} = dsearchn(data.timev_decimated',[10 170]');   % MT7


    % 2) CS-Trials without noise
%             CSIT{1} = [1:5 36:45];                           % MT1
%             CSIT{2} = [1:5 36:45];                           % MT3
%             CSIT{3} = [1 2 3 4 6 8 10 36 39 40 41 42 43 45];  % MT4
%             CSIT{3} = [1:5 36:45]; 
%             CSIT{4} = [1:10 36:45];                           % MT5
%             CSIT{5} = [1:10 36:45];                           % MT6
%             CSIT{6} = [1:10 36:45];                           % MT7

    % 2) CS-Trials without noise for freezing analysis
%              CSIT{1} = [1:5];                           % MT1
%              CSIT{2} = [1:5];                           % MT3
%              CSIT{3} = [1:5];                           % MT4
%              CSIT{4} = [1:5];                           % MT5
%              CSIT{5} = [1:5];                           % MT6
%              CSIT{6} = [1:5];                           % MT7

    %  2) ITI-Trials without noise
%          ITI_ = 10;
%          CSIT_1{1} = 1:10;             % MT1
%          CSIT_1{2} = 1:10;             % MT3
%          CSIT_1{3} = 1:10;             % MT3
%          CSIT_1{3} = [2 3 4 5 8 9] ;   % MT4
%          CSIT_1{4} = 1:10;             % MT5
%          CSIT_1{5} = 1:10;             % MT6
%          CSIT_1{6} = 1:10;             % MT7

    %  2) ITI-Trials without noise
%          ITI_ = 10;
%           CSIT_1{1} = 1:10;             % MT1
%           CSIT_1{2} = 1:10;             % MT3
%           CSIT_1{3} = 1:10;             % MT4
%           CSIT_1{4} = 1:10;             % MT5
%           CSIT_1{5} = 1:10;             % MT6
%           CSIT_1{6} = 1:10;             % MT7


    % Scripts
    %     p_welch_Extinct_full_trials
    %     p_welch_Extinct_full_trials_final_plots
    %     Coherence_Extinct_Retrieval_full_trials
    %     Coherence_Extinct_Retrieval_full_trials_stats
    %     Coherence_Extinct_Retrieval_full_trials_stats_final_plots
    %     Hilbert_phase_Extinct_Full_Trials
    %     Hilbert_phase_Extinct_Full_Trials_final_plots
%         X_corr_Extinct_Full_Trials
%         X_corr_Extinct_full_trials_final_plots
%         Granger_mvgc_full_trials_Extinct_organization
%         Granger_mvgc_full_trials_Extinct_final_plots
%          p_welch_Extinct_Retrieval_behavior
%          p_welch_Extinct_Retrieval_behavior_stats
%          p_welch_Extinct_Retrieval_behavior_stats_final_plots
%         p_welch_Extinct_Retrieval_behavior_fooof
%         p_welch_Extinct_Retrieval_behavior_stats_fooof_final_plots
%         Granger_mvgc_behavior
%         Granger_mvgc_behavior_stats
%         Granger_mvgc_behavior_final_plots

    %% Retrieval Session

    % Load. Choose data from data.lfp according pre_processing.m and main plots script
%        name = strcat(Path,'/',id);
%        fprintf('\n loading data... \n');
% 
%       load(name,'data','parameters');
%      load(name);


    % Baseline
    % This follows the order of animals according to the main loop over ms values (see Analysis.m). Each session needs to be double-checked.

    % - IMPORTANT -
    % The baseline is being analyzed without noise or freezing periods.
    % All CS and ITI are being saved entirely. Only the plots are being saved without noise trials.

    % 1) Baseline timeepochs -> first minute. All animals checked. No noise
%         B_clean{1} = dsearchn(data.timev_decimated',[10 170]');   % MT1
%         B_clean{2} = dsearchn(data.timev_decimated',[10 170]');   % MT3
%         B_clean{3} = dsearchn(data.timev_decimated',[10 170]');   % MT4
%         B_clean{4} = dsearchn(data.timev_decimated',[10 170]');   % MT5
%         B_clean{5} = dsearchn(data.timev_decimated',[10 170]');   % MT6
%         B_clean{6} = dsearchn(data.timev_decimated',[10 170]');   % MT7

    % 1) Baseline timeepochs without noise and freezing behavior
    %       B_clean{1} = dsearchn(data.timev_decimated',[10 60]');    % MT1
    %       B_clean{2} = dsearchn(data.timev_decimated',[60 110]');   % MT3
    %       B_clean{3} = dsearchn(data.timev_decimated',[20 70]');    % MT4
    %       B_clean{4} = dsearchn(data.timev_decimated',[25 75]');    % MT5
    %       B_clean{5} = dsearchn(data.timev_decimated',[10 60]');    % MT6
    %       B_clean{6} = dsearchn(data.timev_decimated',[25 75]');    % MT7

    % 3) Baseline timeepochs -> first minute. All animals checked. No noise
    %       B_clean{1} = dsearchn(data.timev_decimated',[0 60]');    % MT1
    %       B_clean{2} = dsearchn(data.timev_decimated',[0 60]');    % MT3
    %       B_clean{3} = dsearchn(data.timev_decimated',[0 60]');    % MT4
    %       B_clean{4} = dsearchn(data.timev_decimated',[0 60]');    % MT5
    %       B_clean{5} = dsearchn(data.timev_decimated',[0 60]');    % MT6
    %       B_clean{6} = dsearchn(data.timev_decimated',[0 60]');    % MT7

%     % CS-Trials without noise
%           CSIT{1} = 1:5;      % MT1
%           CSIT{2} = 1:5;      % MT3
%           CSIT{3} = 1:5;      % MT4
%           CSIT{3} = 1:5;      % MT4
%           CSIT{4} = 1:5;      % MT5
%           CSIT{5} = 1:5;      % MT6
%           CSIT{6} = 1:5;      % MT7

    % ITI-Trials without noise
    %       ITI_ = 5;
%            CSIT_1{1} = 1:5;      % MT1
%            CSIT_1{2} = 1:5;  % MT3
%           CSIT_1{2} = 1:5;  % MT3
%            CSIT_1{3} = 1:5;    % MT4
%            CSIT_1{3} = 1:5;      % MT4
%            CSIT_1{4} = 1:5;      % MT5
%            CSIT_1{5} = 1:5;      % MT6
%            CSIT_1{6} = 1:5;      % MT7


    % Scripts
%          p_welch_Retrieval_full_trials
%          p_welch_Retrieval_full_trials_final_plots
    %     Hilbert_phase_Retrieval_Full_Trials
    %     Hilbert_phase_Retrieval_Full_Trials_final_plots
%         X_corr_Retrieval_Full_Trials
%         X_corr_Retrieval_full_trials_final_plots
%         Granger_mvgc_full_trials_Retrieval_organization
%         Granger_mvgc_full_trials_Retrieval_final_plots
%         p_welch_Extinct_Retrieval_behavior
%         p_welch_Extinct_Retrieval_behavior_stats
%         p_welch_Extinct_Retrieval_behavior_fooof
%         p_welch_Extinct_Retrieval_behavior_stats_fooof_final_plots
%         p_welch_Extinct_Retrieval_behavior_stats_final_plots
%         Granger_mvgc_behavior
%         Granger_mvgc_behavior_stats
%         Granger_mvgc_behavior_final_plots

    %% Clear
    if ms < length(files.FilesLoaded{1, 1})
        %          clear('FilesLoaded','Path','data','parameters','newStr1','path', 'name' )
        clear('FilesLoaded','Path','newStr1','path', 'name' )

    else
        clear('id','FilesLoaded','Path','ms')
    end

end
toc

fprintf('\n Done. \n');

%% last update 06/09/2024
%  listening:
