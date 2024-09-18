%%  Magnitude Square Coherence  - Organize data to plot and Stats - all animals

% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  03/2024
% Last update: 03/2024


%% Extract data to plot and to stats

% data organization from coher_.mscohr_behavior.F_spect --> see data_preprocessing.m
% - Row 1: Baseline freezing
% - Row 2: Baseline not freezing
% - Row 3: CS-TONE freezing
% - Row 4: CS-TONE not freezing
% - Row 5: ITI freezing
% - Row 6: ITI not freezing


steps                                  = diff(coher_.mscohr_behavior.params.freq); % according to the fft time window

coher_.mscohr_behavior.params.frex_2_12Hz     = 2:steps(1):12;
coher_.mscohr_behavior.params.frex_idx_2_12Hz = dsearchn(coher_.mscohr_behavior.params.freq,coher_.mscohr_behavior.params.frex_2_12Hz');

coher_.mscohr_behavior.params.frex_2_4Hz     = 2:steps(1):4;
coher_.mscohr_behavior.params.frex_idx_2_4Hz = dsearchn(coher_.mscohr_behavior.params.freq,coher_.mscohr_behavior.params.frex_2_4Hz');

coher_.mscohr_behavior.params.frex_2_5_4_5Hz     = 2.5:steps(1):4.5;
coher_.mscohr_behavior.params.frex_idx_2_5_4_5Hz = dsearchn(coher_.mscohr_behavior.params.freq,coher_.mscohr_behavior.params.frex_2_5_4_5Hz');

coher_.mscohr_behavior.params.frex_4_6Hz     = 4:steps(1):6;
coher_.mscohr_behavior.params.frex_idx_4_6Hz = dsearchn(coher_.mscohr_behavior.params.freq,coher_.mscohr_behavior.params.frex_4_6Hz');

coher_.mscohr_behavior.params.frex_6_8Hz     = 6:steps(1):8;
coher_.mscohr_behavior.params.frex_idx_6_8Hz = dsearchn(coher_.mscohr_behavior.params.freq,coher_.mscohr_behavior.params.frex_6_8Hz');

coher_.mscohr_behavior.params.frex_8_10Hz     = 8:steps(1):10;
coher_.mscohr_behavior.params.frex_idx_8_10Hz = dsearchn(coher_.mscohr_behavior.params.freq,coher_.mscohr_behavior.params.frex_8_10Hz');


coher_.mscohr_behavior.stats.full_spectrum_mean   = [];
coher_.mscohr_behavior.stats.spectrum_2_12Hz_mean = [];

coher_.mscohr_behavior.stats.spectrum_2_4Hz_mean      = [];
coher_.mscohr_behavior.stats.spectrum_2_5_4_5Hz_mean  = [];
coher_.mscohr_behavior.stats.spectrum_4_6Hz_mean  = [];
coher_.mscohr_behavior.stats.spectrum_6_8Hz_mean  = [];
coher_.mscohr_behavior.stats.spectrum_8_10Hz_mean  = [];


coher_.mscohr_behavior.stats.spectrum_2_4Hz_peak  = [];
coher_.mscohr_behavior.stats.spectrum_2_5_4_5Hz_peak  = [];
coher_.mscohr_behavior.stats.spectrum_4_6Hz_peak  = [];
coher_.mscohr_behavior.stats.spectrum_6_8Hz_peak  = [];
coher_.mscohr_behavior.stats.spectrum_8_10Hz_peak  = [];


for ii = 1:size(coher_.mscohr_behavior.data,1)
    
    if ii == 1 || ii == 2
        coher_.mscohr_behavior.stats.full_spectrum_mean{ii,1}    = mean(cat(3,coher_.mscohr_behavior.data{ii,1}),3); % baseline

    elseif ii == 3 || ii == 4
        coher_.mscohr_behavior.stats.full_spectrum_mean{ii,1}    = mean(cat(3,coher_.mscohr_behavior.data{ii,CSIT{ms}}),3); % CS-Tone
    else
        coher_.mscohr_behavior.stats.full_spectrum_mean{ii,1}    = mean(cat(3,coher_.mscohr_behavior.data{ii,CSIT_1{ms}}),3); % ITI

    end


    if isempty(coher_.mscohr_behavior.stats.full_spectrum_mean{ii,1})
        coher_.mscohr_behavior.stats.full_spectrum_mean{ii,1} = NaN(3,length(coher_.mscohr_behavior.params.freq));
    end
    

    %coher_.mscohr_behavior.stats.spectrum_2_12Hz_mean{ii,1}  = coher_.mscohr_behavior.stats.full_spectrum_mean{ii,1}(:,coher_.mscohr_behavior.params.frex_idx_2_12Hz);
    coher_.mscohr_behavior.stats.spectrum_2_12Hz_mean{ii,2}  = rescale(coher_.mscohr_behavior.stats.full_spectrum_mean{ii,1}(:,coher_.mscohr_behavior.params.frex_idx_2_12Hz));
   
    coher_.mscohr_behavior.stats.spectrum_2_4Hz_mean{ii,2}       = mean(rescale(coher_.mscohr_behavior.stats.full_spectrum_mean{ii,1}(:,coher_.mscohr_behavior.params.frex_idx_2_4Hz)),2);
    coher_.mscohr_behavior.stats.spectrum_2_5_4_5Hz_mean{ii,2}   = mean(rescale(coher_.mscohr_behavior.stats.full_spectrum_mean{ii,1}(:,coher_.mscohr_behavior.params.frex_idx_2_5_4_5Hz)),2);    
    coher_.mscohr_behavior.stats.spectrum_4_6Hz_mean{ii,2}       = mean(rescale(coher_.mscohr_behavior.stats.full_spectrum_mean{ii,1}(:,coher_.mscohr_behavior.params.frex_idx_4_6Hz)),2);
    coher_.mscohr_behavior.stats.spectrum_6_8Hz_mean{ii,2}       = mean(rescale(coher_.mscohr_behavior.stats.full_spectrum_mean{ii,1}(:,coher_.mscohr_behavior.params.frex_idx_6_8Hz)),2);
    coher_.mscohr_behavior.stats.spectrum_8_10Hz_mean{ii,2}      = mean(rescale(coher_.mscohr_behavior.stats.full_spectrum_mean{ii,1}(:,coher_.mscohr_behavior.params.frex_idx_8_10Hz)),2);

    coher_.mscohr_behavior.stats.spectrum_2_4Hz_peak{ii,2}       = max(rescale(coher_.mscohr_behavior.stats.full_spectrum_mean{ii,1}(:,coher_.mscohr_behavior.params.frex_idx_2_4Hz)),[],2);
    coher_.mscohr_behavior.stats.spectrum_2_5_4_5Hz_peak{ii,2}   = max(rescale(coher_.mscohr_behavior.stats.full_spectrum_mean{ii,1}(:,coher_.mscohr_behavior.params.frex_idx_2_5_4_5Hz)),[],2);
    coher_.mscohr_behavior.stats.spectrum_4_6Hz_peak{ii,2}       = max(rescale(coher_.mscohr_behavior.stats.full_spectrum_mean{ii,1}(:,coher_.mscohr_behavior.params.frex_idx_4_6Hz)),[],2);
    coher_.mscohr_behavior.stats.spectrum_6_8Hz_peak{ii,2}       = max(rescale(coher_.mscohr_behavior.stats.full_spectrum_mean{ii,1}(:,coher_.mscohr_behavior.params.frex_idx_6_8Hz)),[],2);
    coher_.mscohr_behavior.stats.spectrum_8_10Hz_peak{ii,2}      = max(rescale(coher_.mscohr_behavior.stats.full_spectrum_mean{ii,1}(:,coher_.mscohr_behavior.params.frex_idx_8_10Hz)),[],2);

end

clear('ii','steps ')



%% Extract data for each event

coher_.mscohr_behavior.stats.full_spectrum_event           = [];
coher_.mscohr_behavior.stats.spectrum_2_12Hz_event         = [];

coher_.mscohr_behavior.stats.spectrum_2_4Hz_event          = [];
coher_.mscohr_behavior.stats.spectrum_2_5_4_5Hz_event      = [];
coher_.mscohr_behavior.stats.spectrum_4_6Hz_event          = [];
coher_.mscohr_behavior.stats.spectrum_6_8Hz_event          = [];
coher_.mscohr_behavior.stats.spectrum_8_10Hz_event         = [];


coher_.mscohr_behavior.stats.spectrum_2_4Hz_peak_event       = [];
coher_.mscohr_behavior.stats.spectrum_2_5_4_5Hz_peak_event   = [];
coher_.mscohr_behavior.stats.spectrum_4_6Hz_peak_event       = [];
coher_.mscohr_behavior.stats.spectrum_6_8Hz_peak_event       = [];
coher_.mscohr_behavior.stats.spectrum_8_10Hz_peak_event      = [];


for ii = 1:size(coher_.mscohr_behavior.data,1)
    for jj = 1:size(coher_.mscohr_behavior.data,2)


        if isempty(coher_.mscohr_behavior.data{ii,jj})
            continue
        end

        for ll = 1:size(coher_.mscohr_behavior.data{ii,jj},3)
            coher_.mscohr_behavior.stats.spectrum_2_12Hz_event{ii,jj}  = rescale(coher_.mscohr_behavior.data{ii,jj}(:,coher_.mscohr_behavior.params.frex_idx_2_12Hz));

            coher_.mscohr_behavior.stats.spectrum_2_4Hz_event{ii,jj}(:,ll)       = squeeze(mean(rescale(coher_.mscohr_behavior.data{ii,jj}(:,coher_.mscohr_behavior.params.frex_idx_2_4Hz,ll)),2));
            coher_.mscohr_behavior.stats.spectrum_2_5_4_5Hz_event{ii,jj}(:,ll)   = squeeze(mean(rescale(coher_.mscohr_behavior.data{ii,jj}(:,coher_.mscohr_behavior.params.frex_idx_2_5_4_5Hz,ll)),2));
            coher_.mscohr_behavior.stats.spectrum_4_6Hz_event{ii,jj}(:,ll)       = squeeze(mean(rescale(coher_.mscohr_behavior.data{ii,jj}(:,coher_.mscohr_behavior.params.frex_idx_4_6Hz,ll)),2));
            coher_.mscohr_behavior.stats.spectrum_6_8Hz_event{ii,jj}(:,ll)       = squeeze(mean(rescale(coher_.mscohr_behavior.data{ii,jj}(:,coher_.mscohr_behavior.params.frex_idx_6_8Hz,ll)),2));
            coher_.mscohr_behavior.stats.spectrum_8_10Hz_event{ii,jj}(:,ll)      = squeeze(mean(rescale(coher_.mscohr_behavior.data{ii,jj}(:,coher_.mscohr_behavior.params.frex_idx_8_10Hz,ll)),2));

            coher_.mscohr_behavior.stats.spectrum_2_4Hz_peak_event{ii,jj}(:,ll)       = squeeze(max(rescale(coher_.mscohr_behavior.data{ii,jj}(:,coher_.mscohr_behavior.params.frex_idx_2_4Hz,ll)),[],2));
            coher_.mscohr_behavior.stats.spectrum_2_5_4_5Hz_peak_event{ii,jj}(:,ll)   = squeeze(max(rescale(coher_.mscohr_behavior.data{ii,jj}(:,coher_.mscohr_behavior.params.frex_idx_2_5_4_5Hz,ll)),[],2));
            coher_.mscohr_behavior.stats.spectrum_4_6Hz_peak_event{ii,jj}(:,ll)       = squeeze(max(rescale(coher_.mscohr_behavior.data{ii,jj}(:,coher_.mscohr_behavior.params.frex_idx_4_6Hz,ll)),[],2));
            coher_.mscohr_behavior.stats.spectrum_6_8Hz_peak_event{ii,jj}(:,ll)       = squeeze(max(rescale(coher_.mscohr_behavior.data{ii,jj}(:,coher_.mscohr_behavior.params.frex_idx_6_8Hz,ll)),[],2));
            coher_.mscohr_behavior.stats.spectrum_8_10Hz_peak_event{ii,jj}(:,ll)      = squeeze(max(rescale(coher_.mscohr_behavior.data{ii,jj}(:,coher_.mscohr_behavior.params.frex_idx_8_10Hz,ll)),[],2));
        end
    
    end
end

clear('ii','steps ')

%%
figure

% choose data:
not_Frz = 3;
not_NonFrz = 4;

if not_Frz == 1
    sgtitle('Baseline')
elseif not_Frz == 3
    sgtitle('CS-Tone')
else
    sgtitle('ITI')
end



subplot(1,3,1)
if isempty(coher_.mscohr_behavior.stats.spectrum_2_12Hz_mean{not_Frz,2})
   plot(coher_.mscohr_behavior.params.freq(coher_.mscohr_behavior.params.frex_idx_2_12Hz),NaN(1,length(coher_.mscohr_behavior.params.freq(coher_.mscohr_behavior.params.frex_idx_2_12Hz))))
else
   plot(coher_.mscohr_behavior.params.freq(coher_.mscohr_behavior.params.frex_idx_2_12Hz),coher_.mscohr_behavior.stats.spectrum_2_12Hz_mean{not_Frz,2}(1,:))
end

   hold on

if isempty(coher_.mscohr_behavior.stats.spectrum_2_12Hz_mean{not_NonFrz,2})
   plot(coher_.mscohr_behavior.params.freq(coher_.mscohr_behavior.params.frex_idx_2_12Hz),NaN(1,length(coher_.mscohr_behavior.params.freq(coher_.mscohr_behavior.params.frex_idx_2_12Hz))))
else
   plot(coher_.mscohr_behavior.params.freq(coher_.mscohr_behavior.params.frex_idx_2_12Hz),coher_.mscohr_behavior.stats.spectrum_2_12Hz_mean{not_NonFrz,2}(1,:))
end 

ylim([0 1])
legend('Freezing', 'Non-Freezing')


subplot(1,3,2)
if isempty(coher_.mscohr_behavior.stats.spectrum_2_12Hz_mean{not_Frz,2})
   plot(coher_.mscohr_behavior.params.freq(coher_.mscohr_behavior.params.frex_idx_2_12Hz),NaN(1,length(coher_.mscohr_behavior.params.freq(coher_.mscohr_behavior.params.frex_idx_2_12Hz))))
else
   plot(coher_.mscohr_behavior.params.freq(coher_.mscohr_behavior.params.frex_idx_2_12Hz),coher_.mscohr_behavior.stats.spectrum_2_12Hz_mean{not_Frz,2}(2,:))
end

   hold on

if isempty(coher_.mscohr_behavior.stats.spectrum_2_12Hz_mean{not_NonFrz,2})
   plot(coher_.mscohr_behavior.params.freq(coher_.mscohr_behavior.params.frex_idx_2_12Hz),NaN(1,length(coher_.mscohr_behavior.params.freq(coher_.mscohr_behavior.params.frex_idx_2_12Hz))))
else
   plot(coher_.mscohr_behavior.params.freq(coher_.mscohr_behavior.params.frex_idx_2_12Hz),coher_.mscohr_behavior.stats.spectrum_2_12Hz_mean{not_NonFrz,2}(2,:))
end

ylim([0 1])
legend('Freezing', 'Non-Freezing')


subplot(1,3,3)
if isempty(coher_.mscohr_behavior.stats.spectrum_2_12Hz_mean{not_Frz,2})
   plot(coher_.mscohr_behavior.params.freq(coher_.mscohr_behavior.params.frex_idx_2_12Hz),NaN(1,length(coher_.mscohr_behavior.params.freq(coher_.mscohr_behavior.params.frex_idx_2_12Hz))))
else
   plot(coher_.mscohr_behavior.params.freq(coher_.mscohr_behavior.params.frex_idx_2_12Hz),coher_.mscohr_behavior.stats.spectrum_2_12Hz_mean{not_Frz,2}(3,:))
end

   hold on

if isempty(coher_.mscohr_behavior.stats.spectrum_2_12Hz_mean{not_NonFrz,2})
   plot(coher_.mscohr_behavior.params.freq(coher_.mscohr_behavior.params.frex_idx_2_12Hz),NaN(1,length(coher_.mscohr_behavior.params.freq(coher_.mscohr_behavior.params.frex_idx_2_12Hz))))
else
   plot(coher_.mscohr_behavior.params.freq(coher_.mscohr_behavior.params.frex_idx_2_12Hz),coher_.mscohr_behavior.stats.spectrum_2_12Hz_mean{not_NonFrz,2}(3,:))
end

ylim([0 1])
legend('Freezing', 'Non-Freezing')



%% Save
%newStr = regexprep(files.id.name,'.mat','_');
%newStr = files.id(ms).name(1:end-8);
newStr = id(1:end-8);

path = '/Users/flavio/Desktop';
%path = files.FilesLoaded{1, 1}.folder;

name = strcat(path,'/',newStr,'_coher_mscohrwin_1000_nFFT_32768_behavior_1secondTimeW');

% save figure
save(name,'coher_','-v7.3')

% save figure
saveas(gcf,name,'png')

close all

clear('name','newStr','path') 


%% last update 30/03/2024 - 14:01
% listening: