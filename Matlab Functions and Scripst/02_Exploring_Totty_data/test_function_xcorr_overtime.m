function [xcorr_result_full,xcorr_result_full_lag_peak,xcorr_result,xcorr_result_lag_peaks,] = x_corr_overlap(signal1, signal2, nlags, Fs, window_length, overlap, method,fp);

% Cross correlation function over time (hann windowing)
% For non-stationary signals, was considered a window to analyze segments of the signal.

% Important: default parameters have not been programmed yet

% inputs
% signal1 = LFP1
% signal2 = LFP2
% nlags   = number of lags. The maximum lag should be a multiple of the signal periods to capture the cyclic nature of the signal within the
%                           frequency band. Since the signal contains components several components from the frequency band choosed, choosing lags
%                           that cover multiple periods of the lower frequency (longer period) will ensure comprehensive analysis.
% Fs      = samplig frequency
% window_length = time epoch in seconds
% overlap =
% method = according to the parameters of the MATLAB function


% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  04/2024
% Last update: 06/2024

%%
% Default Parameters
% nlags =  501;
% Fs = 1000;
% window_length = 0.865;
% overlap = 0;
% method = 'coeff';
% fp = 1;

% Derived parameters
N = size(x_corr.amplitude{1,2}(2,:),2); % Total number of samples
window_samples = window_length * Fs; % Number of samples per window
overlap_samples = window_samples * (overlap./100); % Number of overlapping samples
step_size = window_samples - overlap_samples; % Step size for moving window

% Generate synthetic signals (examples: sine waves with noise)
signal1 = x_corr.amplitude{1,1}(1,:); % Example signal 1
signal2 = x_corr.amplitude{1,1}(2,:); % Example signal 2

% Initialize cross-correlation result
xcorr_result = [];
xcorr_result_lag_peaks = [];

% Loop over the signals with the defined window and overlap
for start_idx = 1:step_size:N-window_samples
    % Define the window
    end_idx = start_idx + window_samples - 1;
    if end_idx > N
        break; % Ensure we don't exceed the signal length
    end

    % Extract the windowed signal segments
    segment1 = signal1(start_idx:end_idx)';
    segment2 = signal2(start_idx:end_idx)';

    h1 = hann(length(segment1),'periodic');
    h2 = hann(length(segment2),'periodic');

    % Perform cross-correlation
    [xcorr_val, xcorr_lags] = xcorr(segment1.*h1, segment2.*h2, nlags, method);
    %[xcorr_val, xcorr_lags] = xcov(tiedrank(segment1.*h1), tiedrank(segment2.*h2), nlags,method);

    xcorr_val = abs(xcorr_val);
    % Store the result (using xcorr_lags as index for demonstration)
    xcorr_result = [xcorr_result; xcorr_val'];
    % xcorr_result = mapminmax(xcorr_result,0,1);
    % convert lags to ms
    xcorr_lags = xcorr_lags * 1000/Fs;

    % identifies index where the crosscorrelation peaks
    [~, xcorr_lag_peak_idx] = max(abs(xcorr_val));
    % identifies the lag at which the crosscorrelation peaks
    xcorr_lag_peak  = xcorr_lags(xcorr_lag_peak_idx);
    % Store the result (using xcorr_lags as index for demonstration)
    xcorr_result_lag_peaks  = [xcorr_result_lag_peaks;xcorr_lag_peak];

end

% Full session
%[xcorr_result_full, xcorr_full_lags] = xcorr(signal1.*h1, signal2.*h2, nlags, method);
xcorr_result_full = mean(xcorr_result,1);
[~,xcorr_result_full_lag_peak_idx] = max(xcorr_result_full);
xcorr_result_full_lag_peak  = xcorr_lags(xcorr_result_full_lag_peak_idx);

%xcorr_result_full = normalize(xcorr_result_full,'range');


%% Plot the results

if fp == 1

    figure;
    set(gcf,'color','w');
    sc = [1,1,500,1200];
    set(gcf, 'Position', sc);

    sgtitle({'Cross-Correlation Over Time ';['lag = ', num2str(xcorr_result_full_lag_peak), ' ms'];[]},'FontSize',14);

    subplot(5,1,1)
    plot(xcorr_lags,xcorr_result_full,'LineWidth',2,'color',[.6 .6 .6])
    hold on
    plot(xcorr_lags(xcorr_result_full_lag_peak_idx),xcorr_result_full(xcorr_result_full_lag_peak_idx),'ro','MarkerEdgeColor','r','MarkerFaceColor','r')
    title({'mPFC lead <-----> dHPC lead',[]})
    ylabel([{'Normalized'};{'Crosscorrelation'}]);
    xlim([-nlags nlags])
    ylim([min(xcorr_result_full) max(xcorr_result_full)+.1])
    xline(0,'k--','LineWidth', 2)
    box off
    legend('','peak lag','FontSize',10)
    legend('boxoff')

    subplot(5,1,2)
    histogram(xcorr_result_lag_peaks','BinWidth',20,'FaceColor',[.6 .6 .6],'Normalization','pdf')
    ylabel([{'PDF'};{'Peak distribuition'}]);
    % a = fitdist(xcorr_result_lag_peak','normal');
    % hold on
    % plot(a)

    xlim([-nlags nlags])
    box off
    %ylim([0 10])
    xline(0,'k--','LineWidth', 2)

    subplot(5,1,[3 5])
    imagesc(xcorr_lags, (1:size(xcorr_result,1))*step_size/Fs, xcorr_result);
    xlim([-nlags nlags])
    set(gca, 'YDir','normal')
    xline(0,'w--','LineWidth', 2)
    c = colorbar('Ticks',[0.01,0.99],'TickLabels',{'min','max'},'location','southoutside');
    %c.Label.String = 'Normalized peak';
    colormap('parula')
    % Color map matplotlib
    % py_path = "~/anaconda3/envs/Python_3_10/bin/python";
    %     Py_map = getPyPlot_cMap('parula', [], [], py_path);
    %     colormap(Py_map)


    xlabel('Lag (ms)');
    ylabel('Time (seconds)');

    yyaxis right
    r=plot(xcorr_result_lag_peaks,(1:size(xcorr_result,1))*step_size/Fs,'wo','MarkerEdgeColor','k','MarkerFaceColor','w');
    xlim([-nlags nlags])
    a = gca; % Get axis
    a.YColor = 'w';
    a.YTick = [];
    legend(r,'peaks','location','southoutside','FontSize',14)
    legend('boxoff')

else
    fp = 0;
end

end
