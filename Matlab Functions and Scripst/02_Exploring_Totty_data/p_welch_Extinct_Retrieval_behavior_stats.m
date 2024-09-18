%%  Welch power spectral density  - Organize data to plot and Stats - all animals

% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  12/2023
% Last update: 02/2024


%% Extract data to plot ans to stats

% data organization from pw.behavior.F_spect --> see data_preprocessing.m
% - Row 1: Baseline freezing
% - Row 2: Baseline not freezing
% - Row 3: CS-TONE freezing
% - Row 4: CS-TONE not freezing
% - Row 5: ITI freezing
% - Row 6: ITI not freezing


steps                                  = diff(pw.behavior.parameters.freq_); % according to the fft time window

pw.behavior.parameters.frex_2_12Hz     = 2:steps(1):12;
pw.behavior.parameters.frex_idx_2_12Hz = dsearchn(pw.behavior.parameters.freq_,pw.behavior.parameters.frex_2_12Hz');

pw.behavior.parameters.frex_3_6Hz     = 3:steps(1):6;
pw.behavior.parameters.frex_idx_3_6Hz = dsearchn(pw.behavior.parameters.freq_,pw.behavior.parameters.frex_3_6Hz');

pw.behavior.parameters.frex_6_8Hz     = 6:steps(1):8;
pw.behavior.parameters.frex_idx_6_8Hz = dsearchn(pw.behavior.parameters.freq_,pw.behavior.parameters.frex_6_8Hz');

pw.behavior.parameters.frex_6_9Hz     = 6:steps(1):9;
pw.behavior.parameters.frex_idx_6_9Hz = dsearchn(pw.behavior.parameters.freq_,pw.behavior.parameters.frex_6_9Hz');


pw.behavior.stats.full_spectrum_mean   = [];
pw.behavior.stats.spectrum_2_12Hz_mean = [];

pw.behavior.stats.spectrum_3_6Hz_mean  = [];
pw.behavior.stats.spectrum_6_8Hz_mean  = [];
pw.behavior.stats.spectrum_6_9Hz_mean  = [];


pw.behavior.stats.spectrum_3_6Hz_peak  = [];
pw.behavior.stats.spectrum_6_8Hz_peak  = [];
pw.behavior.stats.spectrum_6_9Hz_peak  = [];


for ii = 1:size(pw.behavior.Pxx,1)


    if ii == 1 || ii == 2
        pw.behavior.stats.full_spectrum_mean{ii,1}    = mean(cat(3,pw.behavior.Pxx{ii,1}),3,'omitnan'); % baseline

    elseif ii == 3 || ii == 4
        pw.behavior.stats.full_spectrum_mean{ii,1}    = mean(cat(3,pw.behavior.Pxx{ii,1:5}),3,'omitnan'); % CS-Tone
        %pw.behavior.stats.full_spectrum_mean{ii,1}    = mean(cat(3,pw.behavior.Pxx{ii,:}),3,'omitnan'); % CS-Tone
        
    else
        pw.behavior.stats.full_spectrum_mean{ii,1}    = mean(cat(3,pw.behavior.Pxx{ii,1:5}),3,'omitnan'); % ITI
        %pw.behavior.stats.full_spectrum_mean{ii,1}    = mean(cat(3,pw.behavior.Pxx{ii,:}),3,'omitnan'); % ITI

    end


    if isempty(pw.behavior.stats.full_spectrum_mean{ii,1})
        pw.behavior.stats.full_spectrum_mean{ii,1} = NaN(3,length(pw.behavior.parameters.freq_));
    end
    

    pw.behavior.stats.spectrum_2_12Hz_mean{ii,1}  = pw.behavior.stats.full_spectrum_mean{ii,1}(:,pw.behavior.parameters.frex_idx_2_12Hz);
    pw.behavior.stats.spectrum_2_12Hz_mean{ii,2}  = pw.behavior.stats.full_spectrum_mean{ii,1}(:,pw.behavior.parameters.frex_idx_2_12Hz)./sum(pw.behavior.stats.full_spectrum_mean{ii,1}(:,pw.behavior.parameters.frex_idx_2_12Hz),2);
   
    pw.behavior.stats.spectrum_3_6Hz_mean{ii,2}   = 10.*(mean(pw.behavior.stats.full_spectrum_mean{ii,1}(:,pw.behavior.parameters.frex_idx_3_6Hz)./sum(pw.behavior.stats.full_spectrum_mean{ii,1}(:,pw.behavior.parameters.frex_idx_2_12Hz),2),2));    
    pw.behavior.stats.spectrum_6_8Hz_mean{ii,2}   = 10.*(mean(pw.behavior.stats.full_spectrum_mean{ii,1}(:,pw.behavior.parameters.frex_idx_6_8Hz)./sum(pw.behavior.stats.full_spectrum_mean{ii,1}(:,pw.behavior.parameters.frex_idx_2_12Hz),2),2));
    pw.behavior.stats.spectrum_6_9Hz_mean{ii,2}   = 10.*(mean(pw.behavior.stats.full_spectrum_mean{ii,1}(:,pw.behavior.parameters.frex_idx_6_9Hz)./sum(pw.behavior.stats.full_spectrum_mean{ii,1}(:,pw.behavior.parameters.frex_idx_2_12Hz),2),2));

    pw.behavior.stats.spectrum_3_6Hz_peak{ii,2}   = 10.*(max(pw.behavior.stats.full_spectrum_mean{ii,1}(:,pw.behavior.parameters.frex_idx_3_6Hz)./sum(pw.behavior.stats.full_spectrum_mean{ii,1}(:,pw.behavior.parameters.frex_idx_2_12Hz),2),[],2));
    pw.behavior.stats.spectrum_6_8Hz_peak{ii,2}   = 10.*(max(pw.behavior.stats.full_spectrum_mean{ii,1}(:,pw.behavior.parameters.frex_idx_6_8Hz)./sum(pw.behavior.stats.full_spectrum_mean{ii,1}(:,pw.behavior.parameters.frex_idx_2_12Hz),2),[],2));
    pw.behavior.stats.spectrum_6_9Hz_peak{ii,2}   = 10.*(max(pw.behavior.stats.full_spectrum_mean{ii,1}(:,pw.behavior.parameters.frex_idx_6_9Hz)./sum(pw.behavior.stats.full_spectrum_mean{ii,1}(:,pw.behavior.parameters.frex_idx_2_12Hz),2),[],2));

end

clear('ii','steps ')



%% Extract data for each event

pw.behavior.stats.full_spectrum_event           = [];
pw.behavior.stats.spectrum_2_12Hz_event         = [];

pw.behavior.stats.spectrum_3_6Hz_event      = [];
pw.behavior.stats.spectrum_6_8Hz_event      = [];
pw.behavior.stats.spectrum_6_9Hz_event      = [];


pw.behavior.stats.spectrum_3_6Hz_peak_event   = [];
pw.behavior.stats.spectrum_6_8Hz_peak_event   = [];
pw.behavior.stats.spectrum_6_9Hz_peak_event   = [];


for ii = 1:size(pw.behavior.Pxx,1)
    for jj = 1:size(pw.behavior.Pxx,2)


        if isempty(pw.behavior.Pxx{ii,jj})
            continue
        end

        for ll = 1:size(pw.behavior.Pxx{ii,jj},3)
            pw.behavior.stats.spectrum_2_12Hz_event{ii,jj}  = pw.behavior.Pxx{ii,jj}(:,pw.behavior.parameters.frex_idx_2_12Hz)./sum(pw.behavior.Pxx{ii,jj}(:,pw.behavior.parameters.frex_idx_2_12Hz),2);

            pw.behavior.stats.spectrum_3_6Hz_event{ii,jj}(:,ll)   = squeeze(10.*(mean(pw.behavior.Pxx{ii,jj}(:,pw.behavior.parameters.frex_idx_3_6Hz,ll)./sum(pw.behavior.Pxx{ii,jj}(:,pw.behavior.parameters.frex_idx_2_12Hz,ll),2),2)));
            pw.behavior.stats.spectrum_6_8Hz_event{ii,jj}(:,ll)       = squeeze(10.*(mean(pw.behavior.Pxx{ii,jj}(:,pw.behavior.parameters.frex_idx_6_8Hz,ll)./sum(pw.behavior.Pxx{ii,jj}(:,pw.behavior.parameters.frex_idx_2_12Hz,ll),2),2)));
            pw.behavior.stats.spectrum_6_9Hz_event{ii,jj}(:,ll)      = squeeze(10.*(mean(pw.behavior.Pxx{ii,jj}(:,pw.behavior.parameters.frex_idx_6_9Hz,ll)./sum(pw.behavior.Pxx{ii,jj}(:,pw.behavior.parameters.frex_idx_2_12Hz,ll),2),2)));

            pw.behavior.stats.spectrum_3_6Hz_peak_event{ii,jj}(:,ll)   = squeeze(10.*(max(pw.behavior.Pxx{ii,jj}(:,pw.behavior.parameters.frex_idx_3_6Hz,ll)./sum(pw.behavior.Pxx{ii,jj}(:,pw.behavior.parameters.frex_idx_2_12Hz,ll),2),[],2)));
            pw.behavior.stats.spectrum_6_8Hz_peak_event{ii,jj}(:,ll)       = squeeze(10.*(max(pw.behavior.Pxx{ii,jj}(:,pw.behavior.parameters.frex_idx_6_8Hz,ll)./sum(pw.behavior.Pxx{ii,jj}(:,pw.behavior.parameters.frex_idx_2_12Hz,ll),2),[],2)));
            pw.behavior.stats.spectrum_6_9Hz_peak_event{ii,jj}(:,ll)      = squeeze(10.*(max(pw.behavior.Pxx{ii,jj}(:,pw.behavior.parameters.frex_idx_6_9Hz,ll)./sum(pw.behavior.Pxx{ii,jj}(:,pw.behavior.parameters.frex_idx_2_12Hz,ll),2),[],2)));
        end
    
    end
end

clear('ii','steps ')

%%
figure
sgtitle('Baseline')

subplot(1,3,1)
if isempty(pw.behavior.stats.spectrum_2_12Hz_mean{3,2})
   plot(pw.behavior.parameters.freq_(pw.behavior.parameters.frex_idx_2_12Hz),NaN(1,length(pw.behavior.parameters.freq_(pw.behavior.parameters.frex_idx_2_12Hz))))
else
   plot(pw.behavior.parameters.freq_(pw.behavior.parameters.frex_idx_2_12Hz),pw.behavior.stats.spectrum_2_12Hz_mean{3,2}(1,:))
end

   hold on

if isempty(pw.behavior.stats.spectrum_2_12Hz_mean{3,2})
   plot(pw.behavior.parameters.freq_(pw.behavior.parameters.frex_idx_2_12Hz),NaN(1,length(pw.behavior.parameters.freq_(pw.behavior.parameters.frex_idx_2_12Hz))))
else
   plot(pw.behavior.parameters.freq_(pw.behavior.parameters.frex_idx_2_12Hz),pw.behavior.stats.spectrum_2_12Hz_mean{4,2}(1,:))
end   




subplot(1,3,2)
if isempty(pw.behavior.stats.spectrum_2_12Hz_mean{3,2})
   plot(pw.behavior.parameters.freq_(pw.behavior.parameters.frex_idx_2_12Hz),NaN(1,length(pw.behavior.parameters.freq_(pw.behavior.parameters.frex_idx_2_12Hz))))
else
   plot(pw.behavior.parameters.freq_(pw.behavior.parameters.frex_idx_2_12Hz),pw.behavior.stats.spectrum_2_12Hz_mean{3,2}(2,:))
end

   hold on

if isempty(pw.behavior.stats.spectrum_2_12Hz_mean{3,2})
   plot(pw.behavior.parameters.freq_(pw.behavior.parameters.frex_idx_2_12Hz),NaN(1,length(pw.behavior.parameters.freq_(pw.behavior.parameters.frex_idx_2_12Hz))))
else
   plot(pw.behavior.parameters.freq_(pw.behavior.parameters.frex_idx_2_12Hz),pw.behavior.stats.spectrum_2_12Hz_mean{4,2}(2,:))
end



subplot(1,3,3)
if isempty(pw.behavior.stats.spectrum_2_12Hz_mean{3,2})
   plot(pw.behavior.parameters.freq_(pw.behavior.parameters.frex_idx_2_12Hz),NaN(1,length(pw.behavior.parameters.freq_(pw.behavior.parameters.frex_idx_2_12Hz))))
else
   plot(pw.behavior.parameters.freq_(pw.behavior.parameters.frex_idx_2_12Hz),pw.behavior.stats.spectrum_2_12Hz_mean{3,2}(3,:))
end

   hold on

if isempty(pw.behavior.stats.spectrum_2_12Hz_mean{3,2})
   plot(pw.behavior.parameters.freq_(pw.behavior.parameters.frex_idx_2_12Hz),NaN(1,length(pw.behavior.parameters.freq_(pw.behavior.parameters.frex_idx_2_12Hz))))
else
   plot(pw.behavior.parameters.freq_(pw.behavior.parameters.frex_idx_2_12Hz),pw.behavior.stats.spectrum_2_12Hz_mean{4,2}(3,:))
end


xlim([2 12])


%% Save
%newStr = regexprep(files.id.name,'.mat','_');
%newStr = files.id(ms).name(1:end-8);
newStr = id(1:end-8);

path = '/Users/flavio/Desktop';
%path = files.FilesLoaded{1, 1}.folder;

name = strcat(path,'/',newStr,'_pw_Total_Power_win_1000_nFFT_32768_behavior_1secondTimeW');

% save figure
save(name,'pw','-v7.3')

% save figure
saveas(gcf,name,'png')

close all

clear('name','newStr','path') 


%% last update 26/01/2024
%  listening: Sonic Youth - Disconnection Notice