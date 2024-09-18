% Final Plots Coherence

% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  03/2024
% Last update: 03/2024

%% Load data /retrieval

% cell 1 -> baseline
% cell 2 -> CS-Tone
% cell 3 -> ITI

% - Row 1: mPFC PL <--> mPFC IL
% - Row 2: mPFC PL <--> dHPC
% - Row 3: mPFC IL <--> dHPC


% Coherogram

% Frequency vector
freq_v_coherogram      = coher_.coherogram_MT.params.freq;
freq_v_coherogram_idx  = coher_.coherogram_MT.params.freq_idx_2_12;

% data_2_plot_coherogram_ = [];
% timev_coherogram = []
% Coherogram
for ii = 1:size(coher_.coherogram_MT.data{1,1},1)
    data_2_plot_coherogram_{ii,1}(:,:,:,ms) = coher_.coherogram_MT.stats.data_2_12_norm{ii,1};
    timev_coherogram{ii,1} = coher_.coherogram_MT.params.timev{ii,1}; %Time vector
end



% msCoherence

% Frequency vector
freq_v_mscohr      = coher_.mscohr.params.freq;
freq_v_mscohr_idx  = coher_.mscohr.params.freq_idx_2_12;

% data_2_plot_mscohr_spectrum_all = [];
% data_2_plot_mscohr_band_all = [];
% data_2_plot_mscohr_peak_all = [];

for ii = 1:size(coher_.mscohr.data,1)


    data_2_plot_mscohr_spectrum_all{ii,1}(:,:,ms) = coher_.mscohr.stats.data_2_12_trials_mean{ii,1};            % Spectrum between 1 - 12

    data_2_plot_mscohr_band_all{ii,1}(:,ms) = coher_.mscohr.stats.data_6_8_trials_mean{ii,1};       % choiced mean range
    data_2_plot_mscohr_peak_all{ii,1}(:,ms) = coher_.mscohr.stats.data_6_8_peak_mean{ii,1};   % choiced peak range

end


%% Mean and SEM

if ms == length(files.FilesLoaded{1, 1})

    data_2_plot_coherogram__mean = [];

    for ii = 1:size(data_2_plot_coherogram_,1)
        data_2_plot_coherogram__mean{ii,1} = mean(data_2_plot_coherogram_{ii,1},4);

    end


    data_2_plot_mscohr_mean_SEM = [];

    for ii = 1:size(data_2_plot_mscohr_spectrum_all,1)
        data_2_plot_mscohr_mean_SEM{ii,1}  = mean(data_2_plot_mscohr_spectrum_all{ii,1},3,'omitnan');
        data_2_plot_mscohr_mean_SEM{ii,2}  = std(data_2_plot_mscohr_spectrum_all{ii,1},[],3,'omitnan')./size(data_2_plot_mscohr_spectrum_all{ii,1},3);

        % SEM shades
        curve1_ext{ii,1} = data_2_plot_mscohr_mean_SEM{ii, 1} + data_2_plot_mscohr_mean_SEM{ii, 2};
        curve2_ext{ii,1} = data_2_plot_mscohr_mean_SEM{ii, 1} - data_2_plot_mscohr_mean_SEM{ii, 2};
        data_2_plot_mscohr_mean_SEM{ii,3} = cat(2,curve1_ext{ii,1},flip(curve2_ext{ii,1},2));
        data_2_plot_mscohr_mean_SEM{ii,4} = [freq_v_mscohr(freq_v_mscohr_idx)' fliplr(freq_v_mscohr(freq_v_mscohr_idx)')];

    end


    data_2_plot_mscohr_band_mean = [];
    data_2_plot_mscohr_peak_mean = [];

    for ii = 1:size(data_2_plot_mscohr_band_all,1)
        data_2_plot_mscohr_band_mean{ii,1} = mean(squeeze(data_2_plot_mscohr_band_all{ii,1}),2);
        data_2_plot_mscohr_peak_mean{ii,1} = mean(squeeze(data_2_plot_mscohr_peak_all{ii,1}),2);
    end


    %% Plot

    % choose data
    %data2plot = (coher_.mscohr.data{1,1}(ii,closestfreq)-min(coher_.mscohr.data{1,1}(ii,closestfreq),[],2))./max(coher_.mscohr.data{1,1}(ii,closestfreq),[],2);
    session_part = {'baseline','CS-Tone','ITI'};

    for jj = 1:length(session_part)

        figure(jj+4)%('WindowState','maximized')
        set(gcf,'color','w');
        sc = [1,1,1200,1200];
        set(gcf, 'Position', sc);

        sgtitle({[session_part{jj} ' - Coherogram by multi-taper estimation and Magnitute Square Coherence using Welch{s overlapped averaged ','(window = ' num2str(coher_.mscohr.params.timewin./1000) 's' ' - ' 'overlap = ' num2str(coher_.mscohr.params.overlap) '%)'];[]})


        % Plot coherogram

        titles = {'mPFC PL <-> mPFC IL','mPFC PL <-> dHPC','mPFC IL <-> dHPC'};

%         for ii = 1:size(data_2_plot_coherogram_,1)
%             subplot(4,3,ii)
             contourf(timev_coherogram{jj,1},freq_v_coherogram(freq_v_coherogram_idx),squeeze(data_2_plot_coherogram__mean{jj,1}(ii,:,:)),80,'linecolor','none');
% 
%             title([titles{ii}]);
%             xlabel('Time(s)','FontSize',12)
% 
%             if ii == 1
%                 ylabel('Hz','FontSize',12)
%             end
% 
%             if ii == size(data_2_plot_coherogram__mean,1)
%                 colorbar
%             end
% 
%             %clim([0 1])
% 
%         end

        % Plot MSCoherence

        for ii = 1:size(coher_.mscohr.data{1, 1},1)

            subplot(4,3,ii+3)
            fill(data_2_plot_mscohr_mean_SEM{jj,4},data_2_plot_mscohr_mean_SEM{jj,3}(ii,:),'k','FaceAlpha',0.5,'EdgeColor','none')
            hold on
            plot(freq_v_mscohr(freq_v_mscohr_idx),data_2_plot_mscohr_mean_SEM{jj,1}(ii,:),'color',[0.4, 0.4, 0.4], 'LineWidth',2);
            xlim([2 12])
            ylim([0 .5])

            xlabel('Hz','FontSize',12)

            if ii == 1
                ylabel('Coherence','FontSize',12)
            end


        end

        % Coherence matrix
        % Array to coherence matrix



%         msc_matrix = triu(ones(3)); % number of channels. In this case considering bipolar derivation
%         imAlpha = msc_matrix';
%         imAlpha(1:4:9) = 0;
% 
%         msc_matrix(msc_matrix==0) = data_2_plot_mscohr_peak_mean{jj,1};
%         msc_circular = (msc_matrix + msc_matrix') - 1;
% 
%         myLabel = {'mPFC PL', 'mPFC IL','dHPC'}';
% 
%         subplot(4,3,[7 11])
%         imagesc(msc_matrix,'AlphaData',imAlpha)
% 
%         box off
%         colorbar('southoutside')
%         colormap jet
% 
%         xticklabels(myLabel)
%         xticks(1:length(myLabel))
%         xtickangle(45)
% 
%         yticklabels(myLabel)
%         yticks(1:length(myLabel))
%         title('Coherence Map - 6-8Hz Peak Frequency','FontSize',12)
%         clim([0 1])
% 
%         hAxes = gca;     %Axis handle
% 
%         hAxes.XRuler.Axle.LineStyle = 'none';
%         hAxes.YRuler.Axle.LineStyle = 'none';
% 
%         ylabels = get(hAxes, 'YTickLabel');
%         ylabels{1} = '';   %needs to exist but make it empty
%         set(hAxes, 'YTickLabel', ylabels);
% 
%         xlabels = get(hAxes, 'XTickLabel');
%         xlabels{3} = '';   %needs to exist but make it empty
%         set(hAxes, 'XTickLabel', xlabels);
% 
% 
%         % Circular Map
%         subplot(4,3,[9 12])
% 
%         %myColorMap = ones(length(m_v),3) .* [0 0.7000 0.9000]; * for original circularGraph function to chance color
%         circularGraph(msc_circular,msc_circular,'Label',myLabel);
% 


    end

    clear('closestfreq','data2plot_mscohr','ii','sc','steps','titles','hAxes','imAlpha','msc_circular','msc_matrix','myLabel','xlabels','ylabels')

end

%% Save
% 
% %newStr = regexprep(files.id.name,'.mat','_');
% %newStr = files.id(ms).name(1:end-8);
% newStr = id(1:end-8);
% 
% path = '/Users/flavio/Desktop';
% %path = files.FilesLoaded{1, 1}.folder;
% 
% name = strcat(path,'/',session_part{not1},'_coherence_extinction');
% 
% % save figure
% set(gcf,'renderer', 'painters');
% exportgraphics(gcf,strcat(name,'.png'),'Resolution',300)
% exportgraphics(gcf,strcat(name,'.eps'),'Resolution',300)
% %saveas(gcf,name,'png')
% 
% close all
% clear('name','newStr','path')

%%