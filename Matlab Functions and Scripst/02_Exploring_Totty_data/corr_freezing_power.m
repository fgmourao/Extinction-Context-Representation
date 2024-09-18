% Get freezing time events to correlate with power ranges

% number of events
%sum(cellfun(@length,data.events_behavior([3 4],1:11)),1);
% time_freezing_CS{ms,:} = cellfun(@(x) diff(x,[],2),data.events_behavior(3,CSIT{ms}),'UniformOutput',0);
% corr_extic.time_freezing_CS{ms,1} = [cell2mat(cellfun(@(x) x(:), time_freezing_CS{ms,:}, 'uni', 0)')]' ./ parameters.original_srate;

% normalize for 5 s time epoch
%time_freezing_CS(time_freezing_CS>5) = 5;

%%

% corr_extic.spect_freezing_2_5_4_5{ms,1} = cell2mat(pw.behavior.stats.spectrum_2_4Hz_peak_event(3,CSIT{ms}));
% corr_extic.spect_freezing_2_4{ms,1}    = cell2mat(pw.behavior.stats.spectrum_2_5_4_5Hz_peak_event(3,CSIT{ms}));
% corr_extic.spect_freezing_8_10{ms,1}   = cell2mat(pw.behavior.stats.spectrum_8_10Hz_peak_event(3,CSIT{ms}));

%%
 
% number of events
%sum(cellfun(@length,data.events_behavior([3 4],1:11)),1);
time_freezing_CS{ms,:} = cellfun(@(x) diff(x,[],2),data.events_behavior(3,CSIT{ms}),'UniformOutput',0);
corr_retrieval.time_freezing_CS{ms,1} = [cell2mat(cellfun(@(x) x(:), time_freezing_CS{ms,:}, 'uni', 0)')]' ./ parameters.original_srate;

% normalize for 5 s time epoch
%time_freezing_CS(time_freezing_CS>5) = 5;

%%

% corr_retrieval.spect_freezing_2_5_4_5{ms,1} = cell2mat(pw.behavior.stats.spectrum_2_4Hz_peak_event(3,CSIT{ms}));
% corr_retrieval.spect_freezing_2_4{ms,1}    = cell2mat(pw.behavior.stats.spectrum_2_5_4_5Hz_peak_event(3,CSIT{ms}));
% corr_retrieval.spect_freezing_8_10{ms,1}   = cell2mat(pw.behavior.stats.spectrum_8_10Hz_peak_event(3,CSIT{ms}));