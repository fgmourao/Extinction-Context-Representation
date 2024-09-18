%% Granger Causality - Organize data to plot and Stats

% 1) MVGC Multivariate Granger Causality MATLAB toolbox
% hosted at http://www.sussex.ac.uk/sackler/mvgc
% Current version is mvgc_v1.3, last updated March 2022

% L. Barnett and A. K. Seth, Granger causality for state-space models, Phys. Rev. E 91(4) Rapid Communication, 2015.
% L. Barnett and A. K. Seth, "The MVGC Multivariate Granger Causality Toolbox: A new approach to Granger-causal inference", J. Neurosci. Methods 223, pp 50-68, 2014.

% ----------------------------------

% 2) Time Spectrum and Time vector.
% VAR model estimation regression mode (LWR) by Mike X Cohen (mikexcohen@gmail.com)

% - The code relies on the following functions : - grangerX.m 
%                                                - BSMART ToolBox (http://www.sahs.uth.tmc.edu/hliang/software)


% ----------------------------------

% Part of this code was adapted from "mvgc_demo_statespace.m"
% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  01/2024
% Last update: 02/2024


%% Extract data to plot ans to stats

% data organization from mvgc.F_spect --> see data_preprocessing.m
% - Row 1: Baseline freezing
% - Row 2: Baseline not freezing
% - Row 3: CS-TONE freezing
% - Row 4: CS-TONE not freezing
% - Row 5: ITI freezing
% - Row 6: ITI not freezing


steps                          = diff(mvgc.parameters.freqs); % according to the fft time window

mvgc.parameters.frex_2_12Hz     = 2:steps(1):12;
mvgc.parameters.frex_idx_2_12Hz = dsearchn(mvgc.parameters.freqs,mvgc.parameters.frex_2_12Hz');

mvgc.parameters.frex_3_6Hz     = 3:steps(1):6;
mvgc.parameters.frex_idx_3_6Hz = dsearchn(mvgc.parameters.freqs,mvgc.parameters.frex_3_6Hz');

mvgc.parameters.frex_6_9Hz     = 6:steps(1):9;
mvgc.parameters.frex_idx_6_9Hz = dsearchn(mvgc.parameters.freqs,mvgc.parameters.frex_6_9Hz');

mvgc.stats.full_spectrum_mean   = [];
mvgc.stats.spectrum_2_12Hz_mean = [];

mvgc.stats.spectrum_3_6Hz_mean  = [];
mvgc.stats.spectrum_3_6Hz_peak  = [];

mvgc.stats.spectrum_6_9Hz_mean  = [];
mvgc.stats.spectrum_6_9Hz_peak  = [];


for ii = 1:size(mvgc.F_spect,1)
    
    if ii == 1 || ii == 2
        mvgc.stats.full_spectrum_mean{ii,1}    = mean(cat(4,mvgc.F_spect{ii,1}),4); % baseline

    elseif ii == 3 || ii == 4
        mvgc.stats.full_spectrum_mean{ii,1}    = mean(cat(4,mvgc.F_spect{ii,1:5}),4); % CS-Tone
    else
        mvgc.stats.full_spectrum_mean{ii,1}    = mean(cat(4,mvgc.F_spect{ii,1:5}),4); % ITI

    end


    if isempty(mvgc.stats.full_spectrum_mean{ii,1})
        mvgc.stats.full_spectrum_mean{ii,1} = NaN(mvgc.parameters.number_channels,mvgc.parameters.number_channels,length(mvgc.parameters.freqs));
    end
    

    mvgc.stats.spectrum_2_12Hz_mean{ii,1}  = mvgc.stats.full_spectrum_mean{ii,1}(:,:,mvgc.parameters.frex_idx_2_12Hz);
    mvgc.stats.spectrum_2_12Hz_mean{ii,2}  = mvgc.stats.full_spectrum_mean{ii,1}(:,:,mvgc.parameters.frex_idx_2_12Hz)./sum(mvgc.stats.full_spectrum_mean{ii,1}(:,:,mvgc.parameters.frex_idx_2_12Hz),3);
   

    mvgc.stats.spectrum_3_6Hz_mean{ii,1}   = mean(mvgc.stats.spectrum_2_12Hz_mean{ii,1}(:,:,mvgc.parameters.frex_idx_3_6Hz),3);
    mvgc.stats.spectrum_3_6Hz_mean{ii,2}   = mean(mvgc.stats.spectrum_2_12Hz_mean{ii,2}(:,:,mvgc.parameters.frex_idx_3_6Hz),3);

    mvgc.stats.spectrum_3_6Hz_peak{ii,1}   = max(mvgc.stats.spectrum_2_12Hz_mean{ii,1}(:,:,mvgc.parameters.frex_idx_3_6Hz),[],3);
    mvgc.stats.spectrum_3_6Hz_peak{ii,2}   = max(mvgc.stats.spectrum_2_12Hz_mean{ii,2}(:,:,mvgc.parameters.frex_idx_3_6Hz),[],3);

    mvgc.stats.spectrum_6_9Hz_mean{ii,1}   = mean(mvgc.stats.spectrum_2_12Hz_mean{ii,1}(:,:,mvgc.parameters.frex_idx_6_9Hz),3);
    mvgc.stats.spectrum_6_9Hz_mean{ii,2}   = mean(mvgc.stats.spectrum_2_12Hz_mean{ii,2}(:,:,mvgc.parameters.frex_idx_6_9Hz),3);

    mvgc.stats.spectrum_6_9Hz_peak{ii,1}   = max(mvgc.stats.spectrum_2_12Hz_mean{ii,1}(:,:,mvgc.parameters.frex_idx_6_9Hz),[],3);
    mvgc.stats.spectrum_6_9Hz_peak{ii,2}   = max(mvgc.stats.spectrum_2_12Hz_mean{ii,2}(:,:,mvgc.parameters.frex_idx_6_9Hz),[],3);


end

clear('ii','steps ')

%% Plot
% Plot spectral causal graph. Version 2.
% Plot pairwise spectral quantities in |P|, a 3-dim numerical matrix with
% first index representing target ("to"), second index source ("from")
% and third index frequency range - typically spectral causalities

% all possible combinations
Combinations_2 = nchoosek(1:mvgc.parameters.number_channels,2);
Combinations_1 = flip(Combinations_2,2); % all possible combinations

% data to Plot
data_2_plot = mvgc.stats.spectrum_2_12Hz_mean;
freq_v      = mvgc.parameters.freqs(mvgc.parameters.frex_idx_2_12Hz);

figure
sc = [1,1,960,1200];
set(gcf, 'Position', sc);
sgtitle('Pairwise-conditional Granger causality - frequency domain.');

for cc = 1:size(Combinations_1,1)

    subplot(6,3,cc)
    plot(freq_v,squeeze(data_2_plot{1,1}(Combinations_1(cc,1),Combinations_1(cc,2),:)),'linew',2)
    hold
    plot(freq_v,squeeze(data_2_plot{1,1}(Combinations_2(cc,1),Combinations_2(cc,2),:)),'linew',2)

    if          cc == 1
        legend('PL --> IL', 'IL --> PL')
        ylabel({'Baseline - Freezing';'Granger causality'},'FontWeight','bold','FontSize',12)

    elseif      cc == 2
        legend('PL --> HPC', 'HPC --> PL')
    else
        legend('IL --> HPC', 'HPC --> IL')
    end

    xlabel('(Hz)','FontSize',12)
    xlim([2 12])
    %ylim([0 4])

end


for cc = 1:size(Combinations_1,1)

    subplot(6,3,cc+3)
    plot(freq_v,10.*squeeze(data_2_plot{2,1}(Combinations_1(cc,1),Combinations_1(cc,2),:)),'linew',2)
    hold
    plot(freq_v,10.*squeeze(data_2_plot{2,1}(Combinations_2(cc,1),Combinations_2(cc,2),:)),'linew',2)

    if          cc == 1
        legend('PL --> IL', 'IL --> PL')
        ylabel({'Baseline - NonFreezing';'Granger causality'},'FontWeight','bold','FontSize',12)
    elseif      cc == 2
        legend('PL --> HPC', 'HPC --> PL')
    else
        legend('IL --> HPC', 'HPC --> IL')
    end

    xlabel('(Hz)','FontSize',12)
    xlim([2 12])
    %ylim([0 4])
end



for cc = 1:size(Combinations_1,1)

    subplot(6,3,cc+6)
    plot(freq_v,10.*squeeze(data_2_plot{3,1}(Combinations_1(cc,1),Combinations_1(cc,2),:)),'linew',2)
    hold
    plot(freq_v,10.*squeeze(data_2_plot{3,1}(Combinations_2(cc,1),Combinations_2(cc,2),:)),'linew',2)

    if          cc == 1
        legend('PL --> IL', 'IL --> PL')
        ylabel({'CS-Tone Freezing';'Granger causality'},'FontWeight','bold','FontSize',12)

    elseif      cc == 2
        legend('PL --> HPC', 'HPC --> PL')
    else
        legend('IL --> HPC', 'HPC --> IL')
    end

    xlabel('(Hz)','FontSize',12)
    xlim([2 12])
    %ylim([0 4])

end

for cc = 1:size(Combinations_1,1)

    subplot(6,3,cc+9)
    plot(freq_v,10.*squeeze(data_2_plot{4,1}(Combinations_1(cc,1),Combinations_1(cc,2),:)),'linew',2)
    hold
    plot(freq_v,10.*squeeze(data_2_plot{4,1}(Combinations_2(cc,1),Combinations_2(cc,2),:)),'linew',2)

    if          cc == 1
        legend('PL --> IL', 'IL --> PL')
        ylabel({'CS-Tone NonFreezing';'Granger causality'},'FontWeight','bold','FontSize',12)

    elseif      cc == 2
        legend('PL --> HPC', 'HPC --> PL')
    else
        legend('IL --> HPC', 'HPC --> IL')
    end

    xlabel('(Hz)','FontSize',12)
    xlim([2 12])
    %ylim([0 4])

end

for cc = 1:size(Combinations_1,1)

    subplot(6,3,cc+12)
    plot(freq_v,10.*squeeze(data_2_plot{5,1}(Combinations_1(cc,1),Combinations_1(cc,2),:)),'linew',2)
    hold
    plot(freq_v,10.*squeeze(data_2_plot{5,1}(Combinations_2(cc,1),Combinations_2(cc,2),:)),'linew',2)

    if          cc == 1
        legend('PL --> IL', 'IL --> PL')
        ylabel({'ITI Freezing';'Granger causality'},'FontWeight','bold','FontSize',12)

    elseif      cc == 2
        legend('PL --> HPC', 'HPC --> PL')
    else
        legend('IL --> HPC', 'HPC --> IL')
    end

    xlabel('(Hz)','FontSize',12)
    xlim([2 12])
    %ylim([0 4])

end

for cc = 1:size(Combinations_1,1)

    subplot(6,3,cc+15)
    plot(freq_v,10.*squeeze(data_2_plot{6,1}(Combinations_1(cc,1),Combinations_1(cc,2),:)),'linew',2)
    hold
    plot(freq_v,10.*squeeze(data_2_plot{6,1}(Combinations_2(cc,1),Combinations_2(cc,2),:)),'linew',2)

    if          cc == 1
        legend('PL --> IL', 'IL --> PL')
        ylabel({'ITI NonFreezing';'Granger causality'},'FontWeight','bold','FontSize',12)

    elseif      cc == 2
        legend('PL --> HPC', 'HPC --> PL')
    else
        legend('IL --> HPC', 'HPC --> IL')
    end

    xlabel('(Hz)','FontSize',12)
    xlim([2 12])
    %ylim([0 4])

end


%% Save
newStr = id(1:end-8);
%Path    = files.FilesLoaded{1,1}(ms).folder;
Path = '/Users/flavio/Desktop';
%name_1  = strcat(Path,'/',newStr1,'_Granger_Spec_Mean_CS_Trials_2');
name_  = strcat(Path,'/',newStr,'_Granger_ALL_Spec_Mean_Freezing_10');

saveas(gcf,name_,'png') % save figure

close all

clear('ii','jj')
clear('Combinations_2','Combinations_1','sub_idx','steps','data2plot_mean','cc','session','name_','newStr','path','data2plot')

%% Save data

% Settings
%ms = 1;
newStr = id(1:end-8);
%Path    = files.FilesLoaded{1,1}(ms).folder;
Path = '/Users/flavio/Desktop';
%name_1  = strcat(Path,'/',newStr1,'_Granger_Spec_Mean_CS_Trials_2');
name  = strcat(Path,'/',newStr,'_Granger_Freezing_NonFreezing_stats');


% Save data
save(name,'mvgc','-v7.3')

clear('name','newStr','path')

