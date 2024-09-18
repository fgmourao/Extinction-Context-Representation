%%  Coherence - Stats

% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  12/2023
% Last update: 01/2024


%% coherogram_MT Stats

% Set frequency range for coherogram_MT
steps = diff(coher_.coherogram_MT.params.freq);
coher_.coherogram_MT.params.freq2plot = 2:steps(1):12;
coher_.coherogram_MT.params.freq_idx_2_12 = dsearchn(coher_.coherogram_MT.params.freq,coher_.coherogram_MT.params.freq2plot');

coher_.coherogram_MT.stats = [];

for ii = 1:size(coher_.coherogram_MT.params.combinations,1)
    for jj = 1:size(coher_.coherogram_MT.data,1)
           
            if jj == 1
                coher_.coherogram_MT.stats.data_2_12_{jj,1}(ii,:,:) = coher_.coherogram_MT.data{jj,1}(ii,coher_.coherogram_MT.params.freq_idx_2_12,:);
                coher_.coherogram_MT.stats.data_2_12_norm{jj,1}(ii,:,:) = coher_.coherogram_MT.data{jj,1}(ii,coher_.coherogram_MT.params.freq_idx_2_12,:)./sum(coher_.coherogram_MT.data{jj,1}(ii,coher_.coherogram_MT.params.freq_idx_2_12,:),2);

            else
                coher_.coherogram_MT.stats.data_2_12_{jj,1}(ii,:,:) = mean(coher_.coherogram_MT.data{jj,1}(ii,coher_.coherogram_MT.params.freq_idx_2_12,:,1:5),4);                
                coher_.coherogram_MT.stats.data_2_12_norm{jj,1}(ii,:,:) = mean(coher_.coherogram_MT.data{jj,1}(ii,coher_.coherogram_MT.params.freq_idx_2_12,:,1:5)./sum(coher_.coherogram_MT.data{jj,1}(ii,coher_.coherogram_MT.params.freq_idx_2_12,:,1:5),2),4);
            end

    end
end

%% MS Coherence stats 

% All possible channels combinations

% -> row 1 --> mPFC PL <--> mPFC IL
% -> row 2 --> mPFC PL <--> dHPC
% -> row 2 --> mPFC IL <--> dHPC


% Define frequencies
steps = diff(coher_.mscohr.params.freq);

coher_.mscohr.params.freq2save = 2:steps(1):12;
coher_.mscohr.params.freq_idx_2_12 = dsearchn(coher_.mscohr.params.freq,coher_.mscohr.params.freq2save');

coher_.mscohr.params.freq2save = 3:steps(1):6;
coher_.mscohr.params.freq_idx_3_6 = dsearchn(coher_.mscohr.params.freq,coher_.mscohr.params.freq2save');

coher_.mscohr.params.freq2save = 6:steps(1):8;
coher_.mscohr.params.freq_idx_6_8 = dsearchn(coher_.mscohr.params.freq,coher_.mscohr.params.freq2save');

coher_.mscohr.params.freq2save = 7:steps(1):10;
coher_.mscohr.params.freq_idx_7_10 = dsearchn(coher_.mscohr.params.freq,coher_.mscohr.params.freq2save');

coher_.mscohr.stats = [];

for ii = 1:size(coher_.mscohr.params.combinations,1)
    for jj = 1:size(coher_.mscohr.data,1)
        for tt = 1:size(coher_.mscohr.data{jj,1},3)


            coher_.mscohr.stats.data_2_12_{jj,1}(ii,:,tt)          = coher_.mscohr.data{jj, 1}(ii,coher_.mscohr.params.freq_idx_2_12,tt);
            coher_.mscohr.stats.data_2_12_norm{jj,1}(ii,:,tt)      = 100.*(coher_.mscohr.data{jj, 1}(ii,coher_.mscohr.params.freq_idx_2_12,tt)./sum(coher_.mscohr.data{jj, 1}(ii,coher_.mscohr.params.freq_idx_2_12,tt),2));
            %coher_.mscohr.stats.data_2_12_norm{jj,1}(ii,:,tt)      = detrend(coher_.mscohr.data{jj, 1}(ii,coher_.mscohr.params.freq_idx_2_12,tt),1);

            coher_.mscohr.stats.data_3_6_trials{jj,1}(ii,tt)       = mean(coher_.mscohr.data{jj, 1}(ii,coher_.mscohr.params.freq_idx_3_6,tt),2);
            coher_.mscohr.stats.data_3_6_peak_trials{jj,1}(ii,tt)  = max(coher_.mscohr.data{jj, 1}(ii,coher_.mscohr.params.freq_idx_3_6,tt),[],2);

            coher_.mscohr.stats.data_6_8_trials{jj,1}(ii,tt)       = mean(coher_.mscohr.data{jj, 1}(ii,coher_.mscohr.params.freq_idx_6_8,tt),2);
            coher_.mscohr.stats.data_6_8_peak_trials{jj,1}(ii,tt)  = max(coher_.mscohr.data{jj, 1}(ii,coher_.mscohr.params.freq_idx_6_8,tt),[],2);

            coher_.mscohr.stats.data_7_10_trials{jj,1}(ii,tt)      = mean(coher_.mscohr.data{jj, 1}(ii,coher_.mscohr.params.freq_idx_7_10,tt),2);
            coher_.mscohr.stats.data_7_10_peak_trials{jj,1}(ii,tt) = max(coher_.mscohr.data{jj, 1}(ii,coher_.mscohr.params.freq_idx_7_10,tt),[],2);
       
        end
    end

end

% mean choosen trials

for ii = 1:size(coher_.mscohr.params.combinations,1)
    for jj = 1:size(coher_.mscohr.data,1)

        if jj == 1

            % The redundancy for baseline period here will make the final plot and statistics easier later.
            coher_.mscohr.stats.data_2_12_trials_mean{jj,1}(ii,:) = coher_.mscohr.stats.data_2_12_{jj,1}(ii,:);

            coher_.mscohr.stats.data_3_6_trials_mean{jj,1}(ii,:)  = coher_.mscohr.stats.data_3_6_trials{jj,1}(ii,:);
            coher_.mscohr.stats.data_3_6_peak_mean{jj,1}(ii,:)    = coher_.mscohr.stats.data_3_6_peak_trials{jj,1}(ii,1);

            coher_.mscohr.stats.data_6_8_trials_mean{jj,1}(ii,:)  = coher_.mscohr.stats.data_6_8_trials{jj,1}(ii,1);
            coher_.mscohr.stats.data_6_8_peak_mean{jj,1}(ii,:)    = coher_.mscohr.stats.data_6_8_peak_trials{jj,1}(ii,1);

            coher_.mscohr.stats.data_7_10_trials_mean{jj,1}(ii,:) = coher_.mscohr.stats.data_7_10_trials{jj,1}(ii,1);
            coher_.mscohr.stats.data_7_10_peak_mean{jj,1}(ii,:)   = coher_.mscohr.stats.data_7_10_peak_trials{jj,1}(ii,1);
       
        else
            
            coher_.mscohr.stats.data_2_12_trials_mean{jj,1}(ii,:) = mean(coher_.mscohr.stats.data_2_12_{jj,1}(ii,:,1:5),3);

            coher_.mscohr.stats.data_3_6_trials_mean{jj,1}(ii,:)  = mean(coher_.mscohr.stats.data_3_6_trials{jj,1}(ii,1:5),2);
            coher_.mscohr.stats.data_3_6_peak_mean{jj,1}(ii,:)    = mean(coher_.mscohr.stats.data_3_6_peak_trials{jj,1}(ii,1:5),2);

            coher_.mscohr.stats.data_6_8_trials_mean{jj,1}(ii,:)  = mean(coher_.mscohr.stats.data_6_8_trials{jj,1}(ii,1:5),2);
            coher_.mscohr.stats.data_6_8_peak_mean{jj,1}(ii,:)    = mean(coher_.mscohr.stats.data_6_8_peak_trials{jj,1}(ii,1:5),2);

            coher_.mscohr.stats.data_7_10_trials_mean{jj,1}(ii,:) = mean(coher_.mscohr.stats.data_7_10_trials{jj,1}(ii,1:5),2);
            coher_.mscohr.stats.data_7_10_peak_mean{jj,1}(ii,:)   = mean(coher_.mscohr.stats.data_7_10_peak_trials{jj,1}(ii,1:5),2);
        
        end       
    end
end

clear('steps','ii','jj','tt','coher_.coherogram.params.freq_idx_2_12','coher_.mscohr.params.freq_idx_2_12','closestfreq_3','closestfreq_4')


%% Plot Coherogram by Magnitute Square Coherence using Welch.s overlapped averaged. Time blocks 

% % choose data
% %data2plot = (coher_.mscohr.data{1,1}(ii,closestfreq)-min(coher_.mscohr.data{1,1}(ii,closestfreq),[],2))./max(coher_.mscohr.data{1,1}(ii,closestfreq),[],2);
% %data2plot_coherog = (coher_.coherogram.data{2,1}(:,:,:,2) - mean(coher_.coherogram.data{1,1},3))./std(coher_.coherogram.data{1,1},[],3);
% data2plot_coherog = coher_.coherogram.data{2,1}(:,:,:,2);
% timev_coherog     = coher_.coherogram.params.timev{2,1};
% 
% data2plot_mscohr  = coher_.mscohr.data{2,1}(:,:,2);
% data2plot_mscohr_stat = coher_.mscohr.stats.data_6_8_peak{2,1}(:,1);
% 
% 
% % Set frequency range for coherogram
% steps = diff(coher_.coherogram.params.freq);
% coher_.coherogram.params.freq2plot = 2:steps(1):20;
% coher_.coherogram.params.freq_idx_2_12 = dsearchn(coher_.coherogram.params.freq,coher_.coherogram.params.freq2plot');
% 
% % Set frequency range for mscoherence
% steps = diff(coher_.mscohr.params.freq);
% coher_.mscohr.params.freq2plot = 2:steps(1):12;
% coher_.mscohr.params.freq_idx_2_12 = dsearchn(coher_.mscohr.params.freq,coher_.mscohr.params.freq2plot');
% 
% 
% figure%('WindowState','maximized')
% set(gcf,'color','w');
% sc = [1,1,1200,1200];
% set(gcf, 'Position', sc);
% 
% sgtitle({['Baseline - Coherogram by multi-taper estimation and Magnitute Square Coherence using Welch’s overlapped averaged ','(window = ' num2str(coher_.mscohr.params.timewin./1000) 's' ' - ' 'overlap = ' num2str(coher_.mscohr.params.overlap) '%)'];[]})
% 
% 
% % Plot coherogram
% 
% titles = {'mPFC PL <-> mPFC IL','mPFC PL <-> dHPC','mPFC IL <-> dHPC'};
% 
% for ii = 1:size(coher_.coherogram.data,1)
%     subplot(4,3,ii)
%     contourf(timev_coherog,coher_.coherogram.params.freq(coher_.coherogram.params.freq_idx_2_12),squeeze(data2plot_coherog(ii,coher_.coherogram.params.freq_idx_2_12,:)),80,'linecolor','none');
%     %axis([-1 2.93])
%     caxis([0 1])
%     
%     title([titles{ii}]);
%     xlabel('Time(s)','FontSize',12)
% 
%     if ii == 1
%         ylabel('Hz','FontSize',12)
%     end
% 
%     if ii == size(coher_.coherogram.data,1)
%         colorbar
%     end
% 
%     
%     
% end
% 
% % Plot MSCoherence
% 
% for ii = 1:size(coher_.mscohr.data{1, 1},1)
% 
%     subplot(4,3,ii+3)
% 
%     plot(coher_.mscohr.params.freq(coher_.mscohr.params.freq_idx_2_12),data2plot_mscohr(ii,coher_.mscohr.params.freq_idx_2_12),'color',[0.4, 0.4, 0.4], 'LineWidth',2);
%     xlim([2 12])
%     ylim([0 1])
% 
%     xlabel('Hz','FontSize',12)
% 
%     if ii == 1
%         ylabel('Coherence','FontSize',12)
%     end
%     
% 
% end
% 
% % Coherence matrix
% % Array to coherence matrix
% msc_matrix = triu(ones(3)); % number of channels. In this case considering bipolar derivation
% imAlpha = msc_matrix';
% imAlpha(1:4:9) = 0;
% 
% msc_matrix(msc_matrix==0) = data2plot_mscohr_stat; 
% msc_circular = (msc_matrix + msc_matrix') - 1;
% 
% myLabel = {'mPFC PL', 'mPFC IL','dHPC'}';
% 
% subplot(4,3,[7 11])
% imagesc(msc_matrix,'AlphaData',imAlpha)
% 
% box off
% colorbar('southoutside')
% colormap jet
% 
% xticklabels(myLabel)
% xticks(1:length(myLabel))
% xtickangle(45)
% 
% yticklabels(myLabel)
% yticks(1:length(myLabel))
% title('Coherence Map - 6-8Hz Peak Frequency','FontSize',12)
% clim([0 1])
% 
% hAxes = gca;     %Axis handle
% 
% hAxes.XRuler.Axle.LineStyle = 'none';  
% hAxes.YRuler.Axle.LineStyle = 'none';
% 
% ylabels = get(hAxes, 'YTickLabel');
% ylabels{1} = '';   %needs to exist but make it empty
% set(hAxes, 'YTickLabel', ylabels);
% 
% xlabels = get(hAxes, 'XTickLabel');
% xlabels{3} = '';   %needs to exist but make it empty
% set(hAxes, 'XTickLabel', xlabels);
% 
% 
% % Circular Map
% subplot(4,3,[9 12])
% 
% %myColorMap = ones(length(m_v),3) .* [0 0.7000 0.9000]; * for original circularGraph function to chance color
% circularGraph(msc_circular,msc_circular,'Label',myLabel);
% 
% 
% clear('closestfreq','data2plot_mscohr','data2plot_mscohr_stat','ii','sc','steps','titles','hAxes','imAlpha','msc_circular','msc_matrix','myLabel','xlabels','ylabels')
% 

%% Plot coherogram by multi-taper estimation and mscoherence


timev_coherog     = coher_.coherogram_MT.params.timev{1,1};

% Choose data to plot
data2plot_coherog = coher_.coherogram_MT.stats.data_2_12_{1,1};

data2plot_mscohr  = coher_.mscohr.stats.data_2_12_{1,1}; % Baseline
%data2plot_mscohr  = mean(coher_.mscohr.stats.data_2_12_mean2{1,1}(:,:,CSIT{ms}),3);
data2plot_mscohr_stat = coher_.mscohr.stats.data_6_8_peak_trials{1,1};


figure%('WindowState','maximized')
set(gcf,'color','w');
sc = [1,1,1200,1200];
set(gcf, 'Position', sc);

sgtitle({['Baseline - coherogram MT by multi-taper estimation and Magnitute Square Coherence using Welch’s overlapped averaged ','(window = ' num2str(coher_.mscohr.params.timewin./1000) 's' ' - ' 'overlap = ' num2str(coher_.mscohr.params.overlap) '%)'];[]})


% Plot coherogram_MT

titles = {'mPFC PL <-> mPFC IL','mPFC PL <-> dHPC','mPFC IL <-> dHPC'};

for ii = 1:size(coher_.coherogram_MT.data,1)
    subplot(4,3,ii)
    contourf(timev_coherog,coher_.coherogram_MT.params.freq(coher_.coherogram_MT.params.freq_idx_2_12),squeeze(data2plot_coherog(ii,:,:)),80,'linecolor','none');
    %axis([-1 2.93])
    %caxis([0 1])
    
    title([titles{ii}]);
    xlabel('Time(s)','FontSize',12)

    if ii == 1
        ylabel('Hz','FontSize',12)
    end

    if ii == size(coher_.coherogram_MT.data,1)
        colorbar
    end

    
    
end

% Plot MSCoherence

for ii = 1:size(coher_.mscohr.data{1, 1},1)

    subplot(4,3,ii+3)

    plot(coher_.mscohr.params.freq(coher_.mscohr.params.freq_idx_2_12),data2plot_mscohr(ii,:),'color',[0.4, 0.4, 0.4], 'LineWidth',2);
    xlim([2 12])
    %ylim([0 1])

    xlabel('Hz','FontSize',12)

    if ii == 1
        ylabel('Coherence','FontSize',12)
    end
    

end

% Coherence matrix
% Array to coherence matrix
msc_matrix = triu(ones(3)); % number of channels. In this case considering bipolar derivation
imAlpha = msc_matrix';
imAlpha(1:4:9) = 0;

msc_matrix(msc_matrix==0) = data2plot_mscohr_stat; 
msc_circular = (msc_matrix + msc_matrix') - 1;

myLabel = {'mPFC PL', 'mPFC IL','dHPC'}';

subplot(4,3,[7 11])
imagesc(msc_matrix,'AlphaData',imAlpha)

box off
colorbar('southoutside')
colormap jet

xticklabels(myLabel)
xticks(1:length(myLabel))
xtickangle(45)

yticklabels(myLabel)
yticks(1:length(myLabel))
title('Coherence Map - 6-8Hz Peak Frequency','FontSize',12)
clim([0 1])

hAxes = gca;     %Axis handle

hAxes.XRuler.Axle.LineStyle = 'none';  
hAxes.YRuler.Axle.LineStyle = 'none';

ylabels = get(hAxes, 'YTickLabel');
ylabels{1} = '';   %needs to exist but make it empty
set(hAxes, 'YTickLabel', ylabels);

xlabels = get(hAxes, 'XTickLabel');
xlabels{3} = '';   %needs to exist but make it empty
set(hAxes, 'XTickLabel', xlabels);


% Circular Map
subplot(4,3,[9 12])

%myColorMap = ones(length(m_v),3) .* [0 0.7000 0.9000]; * for original circularGraph function to chance color
circularGraph(msc_circular,msc_circular,'Label',myLabel);


clear('closestfreq','data2plot_mscohr','data2plot_mscohr_stat','ii','sc','steps','titles','hAxes','imAlpha','msc_circular','msc_matrix','myLabel','xlabels','ylabels')

%% Save

%newStr = regexprep(files.id.name,'.mat','_');
%newStr = files.id(ms).name(1:end-8);
newStr = id(1:end-8);

path = '/Users/flavio/Desktop';
%path = files.FilesLoaded{1, 1}.folder;

name = strcat(path,'/',newStr,'_coherence_nfft_32768_fullTrial_2048TimeW');

% save data
save(name,'coher_','-v7.3')

% save figure
saveas(gcf,name,'png')

close all
clear('name','newStr','path') 

