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
% Last update: 05/2024

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
    %     Load function
    %     [data,parameters] = Extracting_LFPs_and_events(id,FilesLoaded,Path);

    %     name = strcat(Path,'/',id);
    %     fprintf('\n loading data... \n');
    %     load(name,'data','parameters');

    %% Optogenetic Stim Conterbalanced Order

    % Habituation -> to set plots, pre-processing and future stats

    OptoSTIM{1} = [8 4]; % 8Hz -> 4Hz
    OptoSTIM{2} = [8 4]; % 8Hz -> 4Hz
    OptoSTIM{3} = [4 8]; % 4Hz -> 8Hz
    OptoSTIM{4} = [4 8]; % 4Hz -> 8Hz
    OptoSTIM{5} = [4 8]; % 4Hz -> 8Hz
    OptoSTIM{6} = [8 4]; % 8Hz -> 4Hz
    OptoSTIM{7} = [4 8]; % 4Hz -> 8Hz

    %% Scripts pre-processing

        Pre_processing_mPCF_dHPC_Habituation_noCS_4Hz_8Hz_Stim;
    %    Pre_processing_mPCF_dHPC_Extinction_CS_8Hz_stim
    %    Pre_processing_mPCF_dHPC_Retrieval_CS_NoStim
    %    Pre_processing_mPCF_dHPC_Renewal_noCS_8Hz_Stim

    %    CorCov_mPFC_dHPC       % Correlation and Covariance Matrices betwwen channels
    %    Main_Plots_mPCF_dHPC_  % plot to check raw data

    %% Habituation
    %  no CS-Tones. 4Hz / 8Hz opto stimulation

    % Load. Choose data from data.lfp according pre_processing.m and main plots script
         name = strcat(Path,'/',id);
         fprintf('\n loading data... \n');
    %
         load(name);

    % Baseline
    % This follows the order of animals according to the main loop over {ms} values (see Analysis.m). Each session needs was double-checked.

    % - IMPORTANT -
    % The baseline is being analyzed without noise or freezing periods.
    % All CS and ITI are being entirely saved. Only final plots will be saved without noise trials.

    % 1) Baseline timeepochs without noise and freezing behavior. Last ~60s before first trial

        B_clean{1} = dsearchn(data.timev_decimated',[1 data.events{2,1}(1,1)./parameters.decimated_srate - 1]');  % R01
%         B_clean{2} = dsearchn(data.timev_decimated',[1 data.events{2,1}(1,1)./parameters.decimated_srate - 1]');  % R02
%         B_clean{3} = dsearchn(data.timev_decimated',[1 data.events{2,1}(1,1)./parameters.decimated_srate - 1]');  % R03
%         B_clean{4} = dsearchn(data.timev_decimated',[1 data.events{2,1}(1,1)./parameters.decimated_srate - 1]');  % R04
%         B_clean{5} = dsearchn(data.timev_decimated',[1 data.events{2,1}(1,1)./parameters.decimated_srate - 1]');  % R06
%         B_clean{6} = dsearchn(data.timev_decimated',[1 data.events{2,1}(1,1)./parameters.decimated_srate - 1]');  % R07
%         B_clean{7} = dsearchn(data.timev_decimated',[1 data.events{2,1}(1,1)./parameters.decimated_srate - 1]');  % R08


    %     B_clean{1} = dsearchn(data.timev_decimated',[data.events{2,1}(1,1)./parameters.decimated_srate - 21 data.events{2,1}(1,1)./parameters.decimated_srate - 11]');  % R01
    %     B_clean{2} = dsearchn(data.timev_decimated',[data.events{2,1}(1,1)./parameters.decimated_srate - 21 data.events{2,1}(1,1)./parameters.decimated_srate - 11]');  % R02
    %     B_clean{3} = dsearchn(data.timev_decimated',[data.events{2,1}(1,1)./parameters.decimated_srate - 21 data.events{2,1}(1,1)./parameters.decimated_srate - 11]');  % R03
    %     B_clean{4} = dsearchn(data.timev_decimated',[data.events{2,1}(1,1)./parameters.decimated_srate - 21 data.events{2,1}(1,1)./parameters.decimated_srate - 11]');  % R04
    %     B_clean{5} = dsearchn(data.timev_decimated',[data.events{2,1}(1,1)./parameters.decimated_srate - 21 data.events{2,1}(1,1)./parameters.decimated_srate - 11]');  % R06
    %     B_clean{6} = dsearchn(data.timev_decimated',[data.events{2,1}(1,1)./parameters.decimated_srate - 21 data.events{2,1}(1,1)./parameters.decimated_srate - 11]');  % R07
    %     B_clean{7} = dsearchn(data.timev_decimated',[data.events{2,1}(1,1)./parameters.decimated_srate - 21 data.events{2,1}(1,1)./parameters.decimated_srate - 11]');  % R08

    % First Stim Blocks without noise
            CSIT{1} = [1:5;1:5];             % R01
    %         CSIT{2} = [1:5;1:5];             % R02
    %         CSIT{3} = [1 2;1 2];                % R03
    %         CSIT{4} = [1:5;1:5];             % R04
    %         CSIT{5} = [1:5;1:5];             % R06
    %         CSIT{6} = [1 2;1 2];              % R07
    %         CSIT{7} = [1 2 5;1 2 5];            % R08

    % First ITI Blocks without noise
    %     CSIT_1{1} = [3:9];                  % R01
        CSIT_1{2} = [1:4 6:10];             % R02
        CSIT_1{3} = [1];                    % R03
        CSIT_1{4} = [1 3 4 5 6 7 9 10];     % R04
        CSIT_1{5} = [2:8];                  % R06
        CSIT_1{6} = [1 2 4 5 6 8 9];        % R07
        CSIT_1{7} = [1:7 9 10];             % R08


    % Scripts
     p_welch_Habituation_noCS_4Hz_8Hz_stim
%     Hilbert_phase_Full_Trials_Habituation_noCS_4Hz_8Hz_Stim
%     sFFT_spectrogram_Habituation_noCS_4Hz_8Hz_Stim

    %% Extinction
    % CS-Tones + 8Hz opto stimulation

    % Load. Choose data from data.lfp according pre_processing.m and main plots script
    %      name = strcat(Path,'/',id);
    %      fprintf('\n loading data... \n');
    %
    %      load(name);

    % Baseline
    % This follows the order of animals according to the main loop over {ms} values (see Analysis.m). Each session needs was double-checked.

    % - IMPORTANT -
    % The baseline is being analyzed without noise or freezing periods.
    % All CS and ITI are being entirely saved. Only final plots will be saved without noise trials.

    % 1) Baseline timeepochs without noise and freezing behavior. Last ~60s before first trial

    %     B_clean{1} = dsearchn(data.timev_decimated',[120 data.events{2,1}(1,1)./parameters.decimated_srate - 1]');  % R02
    %     B_clean{2} = dsearchn(data.timev_decimated',[120 data.events{2,1}(1,1)./parameters.decimated_srate - 1]');  % R04

    %     B_clean{1} = dsearchn(data.timev_decimated',[data.events{2,1}(1,1)./parameters.decimated_srate - 21 data.events{2,1}(1,1)./parameters.decimated_srate - 11]');  % R02
    %     B_clean{2} = dsearchn(data.timev_decimated',[data.events{2,1}(1,1)./parameters.decimated_srate - 21 data.events{2,1}(1,1)./parameters.decimated_srate - 11]');  % R04
    %     B_clean{3} = dsearchn(data.timev_decimated',[data.events{2,1}(1,1)./parameters.decimated_srate - 21 data.events{2,1}(1,1)./parameters.decimated_srate - 11]');  % R07

    % First Stim Blocks without noise
    %      CSIT{1} = [1:45];             % R02
    %      CSIT{2} = [1:45];             % R04
    %      CSIT{3} = [1:10];             % R07

    % First ITI Blocks without noise
    %     CSIT_1{1} = [1:4 6:10];             % R02
    %     CSIT_1{2} = [1 3 4 5 6 7 9 10];     % R04
    %     CSIT_1{3} = [1 3 4 5 6 7 9 10];     % R07


    % Scripts
    % p_welch_Extinction_CS_8Hz_stim
    % Hilbert_phase_Full_Trials_Extinction_CS_8Hz_Stim

    %% Retrieval
    %  CS-Tones - No Opto Stimulation

    % Load. Choose data from data.lfp according pre_processing.m and main plots script
    %          name = strcat(Path,'/',id);
    %          fprintf('\n loading data... \n');
    %
    %          load(name);

    % Baseline
    % This follows the order of animals according to the main loop over {ms} values (see Analysis.m). Each session needs was double-checked.

    % - IMPORTANT -
    % The baseline is being analyzed without noise or freezing periods.
    % All CS and ITI are being entirely saved. Only final plots will be saved without noise trials.

    % 1) Baseline timeepochs without noise and freezing behavior. Last ~60s before first trial

    %     B_clean{1} = dsearchn(data.timev_decimated',[data.events{2,1}(1,1)./parameters.decimated_srate - 21 data.events{2,1}(1,1)./parameters.decimated_srate - 11]');  % R02
    %     B_clean{2} = dsearchn(data.timev_decimated',[data.events{2,1}(1,1)./parameters.decimated_srate - 21 data.events{2,1}(1,1)./parameters.decimated_srate - 11]');  % R04
    %     B_clean{3} = dsearchn(data.timev_decimated',[data.events{2,1}(1,1)./parameters.decimated_srate - 21 data.events{2,1}(1,1)./parameters.decimated_srate - 11]');  % R07

    % First Stim Blocks without noise
    %     CSIT{1} = [1:5];             % R02
    %     CSIT{2} = [1:5];             % R04
    %     CSIT{3} = [1:5];             % R07

    % First ITI Blocks without noise
    %     CSIT_1{1} = [1:5];     % R02
    %     CSIT_1{2} = [1:5];     % R04
    %     CSIT_1{3} = [1:5]];    % R07


    %     % Scripts
    %     p_welch_Retrieval_CS_NoStim
    %     Hilbert_phase_Full_Trials_Retrieval_CS_NoStim


    %% Renewal
    %  no CS-Tones. 8Hz opto stimulation

    % Load. Choose data from data.lfp according pre_processing.m and main plots script
%     name = strcat(Path,'/',id);
%     fprintf('\n loading data... \n');
% 
%     load(name);

    % Baseline
    % This follows the order of animals according to the main loop over {ms} values (see Analysis.m). Each session needs was double-checked.

    % - IMPORTANT -
    % The baseline is being analyzed without noise or freezing periods.
    % All CS and ITI are being entirely saved. Only final plots will be saved without noise trials.

    % 1) Baseline timeepochs without noise and freezing behavior. Last ~60s before first trial

%     B_clean{1} = dsearchn(data.timev_decimated',[data.events{2,1}(1,1)./parameters.decimated_srate - 31 data.events{2,1}(1,1)./parameters.decimated_srate - 21]');  % R02
%     B_clean{2} = dsearchn(data.timev_decimated',[data.events{2,1}(1,1)./parameters.decimated_srate - 31 data.events{2,1}(1,1)./parameters.decimated_srate - 21]');  % R04
%     B_clean{3} = dsearchn(data.timev_decimated',[data.events{2,1}(1,1)./parameters.decimated_srate - 31 data.events{2,1}(1,1)./parameters.decimated_srate - 21]');  % R07

    % First Stim Blocks without noise
%     CSIT{1} = [1:10];             % R02
%     CSIT{2} = [1:10];             % R04
%     CSIT{3} = [1:10];             % R07

    % First ITI Blocks without noise
    %     CSIT_1{1} = [1:5];     % R02
    %     CSIT_1{2} = [1:5];     % R04
    %     CSIT_1{3} = [1:5]];    % R07


    % Scripts
%     p_welch_Renewal_noCS_8Hz_Stim
%     Hilbert_phase_Full_Trials_Renewal_noCS_8Hz_stim

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

%% last update 20/05/2024
%  listening:
