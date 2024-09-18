%% Granger - Armor model by MickX

% Choose channels
Combinations_ = nchoosek(1:size(mvgc.data,1),2);

% Define band frequency
% 2 - 12 Hertz
steps         = diff(mvgc.parameters.freqs); % according to the fft time window
mvgc.parameters.frex_2_12Hz     = 2:steps(1):12;
mvgc.parameters.frex_idx_2_12Hz = dsearchn(mvgc.parameters.freqs,mvgc.parameters.frex_2_12Hz');

win = 1000; % in seconds.
mvgc.parameters.time_winx  = round(win/(1000/mvgc.parameters.fs)); % in samples

for tt = 1:size(mvgc.data,1)
    for ii = 1:size(mvgc.data{tt,1},3)

        for cc = 1:size(Combinations_,1)

            [mvgc.Time_x2y_Spec{tt,ii}(:,:,cc),mvgc.Time_y2x_Spec{tt,ii}(:,:,cc),mvgc.Time_x2y_FInt{tt,ii}(:,:,cc),mvgc.Time_y2x_FInt{tt,ii}(:,:,cc)] = ... 
                grangerX(mvgc.data{tt,ii}(Combinations_(cc,:),:,ii),mvgc.parameters.fs,mvgc.parameters.freqs(mvgc.parameters.frex_idx_2_12Hz),mvgc.parameters.time_winx,mvgc.parameters.morder);

        end
   end
end

%% PLot Granger Time Spec. -> Need to fix

% data_2_plot_Spec = movmean(x2yS,225,2);
% time_v_spec = linspace(0,10,size(data_2_plot_Spec,2));
% 
% data_2_plot_Behav = decimate(data.behavior{3,4}(1,1:end-1),4);
% 
% for ii = 1:length(data_2_plot_Behav) - mvgc.parameters.time_winx
%     data_2_plot_Behav_m(1,ii) = mean(data_2_plot_Behav(1,ii:ii+mvgc.parameters.time_winx),2);
% end
% 
% time_v_Behav = linspace(0,10,size(data_2_plot_Behav_m,2));
% 
% % data_2_plot_Freez_idx(1,:) = data.behavior{5, 2}(1,:);
% % data_2_plot_Freez_idx(2,:) = data.behavior{5, 2}(1,:) +  data.behavior{5, 2}(2,:)-1;
% 
% figure
% set(gcf,'color','w');
% 
% %sgtitle({'Amplitude Spectrum via short-window FFT';['Hamming window = ' num2str(short_fft.timewin./1000) 's' ' - ' 'overlap = ' num2str(short_fft.overlap) '%']})
% 
% contourf(time_v_spec, mvgc.parameters.freqs(mvgc.parameters.frex_idx),data_2_plot_Spec,80,'linecolor','none');
% xlabel('Time (s)','FontSize',14), ylabel('Frequency (Hz)','FontSize',14)
% clim([0 1.4])
% c = colorbar;
% c.Label.String = 'Granger (A.U.)';
% 
% hold on
% 
% yyaxis right
% plot(time_v_Behav,data_2_plot_Behav_m,'linew',2,'color',[1 1 1])
% %plot([time_v_Behav(data_2_plot_Freez_idx(1,:));time_v_Behav(data_2_plot_Freez_idx(2,:))], [ones(1,length(data_2_plot_Freez_idx(1,:))).*20;ones(1,length(data_2_plot_Freez_idx(2,:))).*20],'k-','linew', 3,'Color',[1, 0, .2])
% 
% 
% 
% 
% 
% 
% clear('ch','steps','win')