
%% Coherence from behavior epochs

% 1)
% Magnitute Square Coherence using Welch’s overlapped averaged
% Matlab based welch function --> mscoher

% 2)
% LFP coherogram by multi-taper estimation
% Adapted from --> bz_MTCoherogram.m
% Copyright (C) 2010-2014 by Michaël Zugaro

% The code relies on the following functions:
% --> bz_MTCoherogram -> https://github.com/buzsakilab/buzcode/blob/master/analysis/SpectralAnalyses/bz_MTCoherogram.m
% --> Chronux analysis software -> http://chronux.org
% --> Freely Moving Animal (FMA) Toolbox -> http://fmatoolbox.sourceforge.net/FMAToolbox


% Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  03/2024
% Last update: 03/2024

%%
fprintf('\n Magnitude-squared coherence estimate... \n');

%%  mscoher

% variable: coher_ :

% First row: baseline
% Second row: CS tones
% Third row: IT
% 
% In each cell:
%     - Rows: combinations
%     - Columns: frequencies


coher_.mscohr_behavior   = []; % mscoher from matlab built function


% All possible channels combinations
coher_.mscohr_behavior.params.combinations  = nchoosek(1:size(data.lfp{5,1},1),2);
% -> row 1 --> mPFC PL <--> mPFC IL
% -> row 2 --> mPFC PL <--> dHPC
% -> row 2 --> mPFC IL <--> dHPC

% Parameters
% Time window
coher_.mscohr_behavior.params.timewin    = 725; % in ms

% Convert time window to points
coher_.mscohr_behavior.params.timewinpnts  = hamming(round(coher_.mscohr_behavior.params.timewin/(1000/parameters.decimated_srate)));

% Number of overlap samples
coher_.mscohr_behavior.params.overlap = 95;
coher_.mscohr_behavior.params.noverlap = floor(coher_.mscohr_behavior.params.overlap*0.01*(round(coher_.mscohr_behavior.params.timewin/(1000/parameters.decimated_srate))));

% nFFT
coher_.mscohr_behavior.params.nFFT = 2^15; %2^nextpow2(round(coher_.mscohr_behavior.params.timewin/(1000/parameters.srate)));

data_2_use = data.lfp_behavior;

coher_.mscohr_behavior.data = cell(size(data.lfp_behavior));

for ii = 1:size(data_2_use,1)
    for jj = 1:size(data_2_use,2)

        if isempty(data_2_use{ii,jj})
            continue
        end

        for tt = 1:size(data_2_use{ii,jj},3)

            for mm = 1:length(coher_.mscohr_behavior.params.combinations)
                temp_1 = data_2_use{ii,jj}(coher_.mscohr_behavior.params.combinations(mm,1),:,tt);
                temp_2 = data_2_use{ii,jj}(coher_.mscohr_behavior.params.combinations(mm,2),:,tt);

%                 if ii == 1
                    [coher_.mscohr_behavior.data{ii,jj}(mm,:,tt),coher_.mscohr_behavior.params.freq] = mscohere(temp_1,temp_2,coher_.mscohr_behavior.params.timewinpnts,coher_.mscohr_behavior.params.overlap,coher_.mscohr_behavior.params.nFFT,parameters.decimated_srate);
%                 else
%                     coher_.mscohr_behavior.data{ii,jj}(mm,:,tt) = mscohere(temp_1,temp_2,coher_.mscohr_behavior.params.timewinpnts,coher_.mscohr_behavior.params.overlap,coher_.mscohr_behavior.params.nFFT,parameters.decimated_srate);

%                 end
            end
        end
    end
end


clear ('ii','jj','cc','tt','mm','data_2_use','temp_1','temp_2')



%% IMPORTANT
% NOTE ------> Coherogram NEEDS TO BE CHANGE based on behavior time windows
% For now just copied and paste from full trials script .

%% Coherogram
%
%    [coherogram,phase,t,f] = MTCoherogram(lfp1,lfp2,<options>)
%
%    lfp1,lfp2      wide-band LFPs (one channel each).
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'srate'       sampling rate (in Hz) (default = from timestamps if
%                   available, otherwise 1250Hz)
%     'range'       frequency range (in Hz) (default = all)
%     'window'      duration (in s) of the time window (default = 5)
%     'overlap'     overlap between successive windows (default = window/2)
%     'step'        step between successive windows (default = window/2)
%     'tapers'      relative resolution and order of the tapers [NW K]
%                   (default = [3 5])
%     'pad'         FFT padding (see help for <a href="matlab:help cohgramc">cohgramc</a>) (default = 0)
%     'show'        plot results (default = 'off')
%     'cutoffs'     cutoff values for color plot (default = [0 1])
%    =========================================================================
%
%  NOTES
%
%    The LFP can be provided either as a time stamped matrix (list of time-voltage
%    pairs), or as a voltage vector - in which case the frequency must be specified.
%
%    The time displacement between successive short time coherences can be supplied
%    either as a 'step' (explicit time difference) or as an 'overlap' (between
%    successive time windows).
%
%  OUTPUT
%
%    coherogram     coherogram magnitude
%    phase          coherogram phase
%    t              time bins
%    f              frequency bins
%
%%
% fprintf('\n Coherogram by multi-taper estimation... \n');
% 
% startup_chronux
% startup_FMA
% 
% clear("chronux_root",'startup_FMA','FMA_root')

%% Parameters for coherogram following instructions above

% First row: baseline
% In each cell:
%     - Rows: combinations
%     - Columns: frequencies
%     - 3th dimention: time

% Second row: CS tones
% Third row: ITI - Inter trials

% In each cell:
%     - Rows: combinations
%     - Columns: frequencies
%     - 3th dimention: time
%     - 4th dimention: trials


% coher_.coherogram = [];
% 
% coher_.coherogram.params.combinations  = nchoosek(1:size(data.lfp{5,1},1),2);
% coher_.coherogram.params.srate = parameters.decimated_srate;
% coher_.coherogram.params.range = [1 coher_.coherogram.params.srate/2];
% coher_.coherogram.params.window = 5.000;
% 
% overlap = .85; % (%)
% 
% coher_.coherogram.params.overlap = overlap * coher_.coherogram.params.window;
% coher_.coherogram.params.step = coher_.coherogram.params.window - coher_.coherogram.params.overlap;
% coher_.coherogram.params.tapers = [3 6];
% coher_.coherogram.params.pad = 2;
% coher_.coherogram.params.show = 'off';
% coher_.coherogram.params.cutoffs = [0 1];

%% LFP coherogram by multi-taper estimation

% All possible channels combinations following mscoherence parameters defined above

% baseline

% not1 = 6;
% 
% for ii = 1:length(coher_.coherogram.params.combinations)
%     temp_1 = data.lfp{not1,1}(coher_.coherogram.params.combinations(ii,1),B_clean{ms}(1):B_clean{ms}(2));
%     temp_2 = data.lfp{not1,1}(coher_.coherogram.params.combinations(ii,2),B_clean{ms}(1):B_clean{ms}(2));
% 
% 
%     [coher_.coherogram.data{1,1}(ii,:,:),coher_.coherogram.phase{1,1}(ii,:,:),coher_.coherogram.timev{1,1},coher_.coherogram.freqv] = ...
%         bz_MTCoherogram(temp_1',temp_2','frequency',coher_.coherogram.params.srate,'range',coher_.coherogram.params.range,'window',coher_.coherogram.params.window,...
%         'overlap',coher_.coherogram.params.overlap,'step',coher_.coherogram.params.step,'tapers',coher_.coherogram.params.tapers,'pad',coher_.coherogram.params.pad,'show',coher_.coherogram.params.show,'cutoffs',coher_.coherogram.params.cutoffs);
% 
% end
% 
% 
% % CS-tones
% not1 = 7;
% 
% for ii = 1:length(coher_.coherogram.params.combinations)
%     for jj = 1:size(data.lfp{not1,1},3)
% 
%         temp_1 = data.lfp{not1,1}(coher_.coherogram.params.combinations(ii,1),:,jj);
%         temp_2 = data.lfp{not1,1}(coher_.coherogram.params.combinations(ii,2),:,jj);
% 
%         [coher_.coherogram.data{2,1}(ii,:,:,jj),coher_.coherogram.phase{2,1}(ii,:,:,jj),coher_.coherogram.timev{2,1},coher_.coherogram.freqv] = ...
%             bz_MTCoherogram(temp_1',temp_2','frequency',coher_.coherogram.params.srate,'range',coher_.coherogram.params.range,'window',coher_.coherogram.params.window,...
%             'overlap',coher_.coherogram.params.overlap,'step',coher_.coherogram.params.step,'tapers',coher_.coherogram.params.tapers,'pad',coher_.coherogram.params.pad,'show',coher_.coherogram.params.show,'cutoffs',coher_.coherogram.params.cutoffs);
% 
%     end
% end
% 
% % ITI - ( Inter trials)
% not1 = 8;
% 
% for ii = 1:length(coher_.coherogram.params.combinations)
%     for jj = 1:size(data.lfp{not1,1},3)
% 
%         temp_1 = data.lfp{not1,1}(coher_.coherogram.params.combinations(ii,1),:,jj);
%         temp_2 = data.lfp{not1,1}(coher_.coherogram.params.combinations(ii,2),:,jj);
% 
%         [coher_.coherogram.data{3,1}(ii,:,:,jj),coher_.coherogram.phase{3,1}(ii,:,:,jj),coher_.coherogram.timev{3,1},coher_.coherogram.freqv] = ...
%             bz_MTCoherogram(temp_1',temp_2','frequency',coher_.coherogram.params.srate,'range',coher_.coherogram.params.range,'window',coher_.coherogram.params.window,...
%             'overlap',coher_.coherogram.params.overlap,'step',coher_.coherogram.params.step,'tapers',coher_.coherogram.params.tapers,'pad',coher_.coherogram.params.pad,'show',coher_.coherogram.params.show,'cutoffs',coher_.coherogram.params.cutoffs);
% 
%     end
% end
% 
% 
% clear ('ii','jj','temp_1','temp_2','not1','overlap')

%% Save

%newStr = regexprep(files.id.name,'.mat','_');
%newStr = files.id(ms).name(1:end-8);
newStr = id(1:end-8);

path = '/Users/flavio/Desktop';
%path = files.FilesLoaded{1, 1}.folder;

name = strcat(path,'/',newStr,'_coher_mscohrwin_1000_nFFT_32768_behavior_1secondTimeW');

% save data
save(name,'coher_','-v7.3')

clear('name','newStr','path') 

%% last update 30/03/2024 - 14:01
% listening:
