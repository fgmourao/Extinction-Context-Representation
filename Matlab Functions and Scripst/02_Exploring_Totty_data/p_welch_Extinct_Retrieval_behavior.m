%%  Welch power spectral density estimate 

% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  12/2023
% Last update: 01/2024


% First extract the data with: Extracting_LFPs_and_events.m
% ... then organize the data with the script: Pre_processing.m

%%
fprintf('\n Welch power spectral density estimate ... \n');

%% pwelch
%  [pxx,f] = pwelch(x,window,noverlap,f,fs)

% - cell          
% - 1st column    - > baseline
% - 2nd column    - > cs events

% in each cell
% - rows          - > Hz
% - columns       - > channels
% 3th dimentions  - > CS-Trials

pw = [];

% Time window
pw.behavior.parameters.timewin    = 1000; % in ms

% Convert time window to points
pw.behavior.parameters.timewinpnts   = hamming(round(pw.behavior.parameters.timewin/(1000/parameters.decimated_srate)));

% nFFT
pw.behavior.parameters.nFFT = 2^15; %4096; %2^nextpow2(pw.behavior.baseline_timewinpnts));

% Number of overlap samples
pw.behavior.parameters.overlap = 90;
pw.behavior.parameters.noverlap = floor(pw.behavior.parameters.overlap*0.01 * pw.behavior.parameters.timewin);


data_2_use = data.lfp_behavior;

pw.behavior.Pxx = cell(size(data.lfp_behavior));

for ii = 1:size(data_2_use,1)
    for jj = 1:size(data_2_use,2)

        if isempty(data_2_use{ii,jj})
            continue
        end

        for cc = 1:size(data_2_use{ii,jj},1)
            for tt = 1:size(data_2_use{ii,jj},3)

                if ii == 1 || ii == 2
                    [pw.behavior.Pxx{ii,jj}(cc,:,tt),pw.behavior.parameters.freq_] = pwelch(data_2_use{ii,jj}(cc,:,tt), pw.behavior.parameters.timewinpnts, pw.behavior.parameters.noverlap, pw.behavior.parameters.nFFT, parameters.decimated_srate);

                else

                    pw.behavior.Pxx{ii,jj}(cc,:,tt) = pwelch(data_2_use{ii,jj}(cc,:,tt), pw.behavior.parameters.timewinpnts, pw.behavior.parameters.noverlap, pw.behavior.parameters.nFFT, parameters.decimated_srate);
                
                end
            end
        end                
    end

end


clear ('ii','jj','cc','tt','data_2_use')

%%
% figure
% plot(pw.behavior.parameters.freq_,pw.behavior.Pxx{1,1}(1,:))
% hold
% plot(pw.behavior.parameters.freq_,pw.behavior.Pxx{2,1}(1,:))
% 
% xlim([2 12])


%% Save

%newStr = regexprep(files.id.name,'.mat','_');
%newStr = files.id(ms).name(1:end-8);
newStr = id(1:end-8);

path = '/Users/flavio/Desktop';
%path = files.FilesLoaded{1, 1}.folder;

name = strcat(path,'/',newStr,'_pw_Total_Power_win_1000_nFFT_32768_behavior_1secondTimeW');

% save data
save(name,'pw','-v7.3')

clear('name','newStr','path') 

%% last update 04/03/2024
%  listening: Sonic Youth - Disconnection Notice
