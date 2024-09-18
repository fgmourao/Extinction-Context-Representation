
%% Cross-correlation based on Instantaneous Amplitude (Hilbert Transform)
%  - Performs analysis considering the trial periods

% Based on :

% - 10.1016/j.jneumeth.2010.06.019

% - 10.1038/nn.4327

% - Analyzing neural Time Series data - Theory and Practice
%   Mike X Cohen
%   Chapter 27, Power based conectivity pg 361 -
%   Following the model proposed by Mike Cohen, the function x_corr could
%   be eplaced by xcov to scale to the correlation coefficient. Since the xcov function does not perform Spearman's correction,
%   the tiedrank function was applied to the data to perform the rank-transform.

% Relies in the following function:
% - x_corr_overtime

% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  04/2024
% Last update: 06/2024

%%
% First extract the data with: Extracting_LFPs_and_events.m
% ... then organize the data with the script: Pre_processing.m

%% Run each session sequentially

%%
fprintf('\n Cross-correlation based on Instantaneous Amplitude... \n');

%% Filtering and Hilbert Transform to extract amplitude

% Loop over channels and make hilbert transform
% Here, periods without noise are considered for both the baseline and the trials.
% Defined by the variables: Baseline: B_clean / CS-Tone/ITI: CSIT /

scale = 1/0.153;

x_corr = [];

x_corr.parameters.filters = [2 5;6 8;6 12];

% Baseline
for ii = 1:size(data.lfp{6,1},1)
    for hh = 1:size(x_corr.parameters.filters,1)

        temp = eegfilt(scale.*(data.lfp{6,1}(ii,B_clean{1,ms}(1,1):B_clean{1,ms}(2,1))),parameters.decimated_srate,x_corr.parameters.filters(hh,1),x_corr.parameters.filters(hh,2)); % just filtering
        x_corr.amplitude{1,hh}(ii,:) = abs(hilbert(temp))-mean(abs(hilbert(temp)),2); % getting the amplitude envelope and subtracted the mean to eliminate each signal's DC component
        % x_corr.amplitude{1,hh}(ii,:) = temp;

    end

end


%  CS-Tone
for ii = 1:size(data.lfp{7,1},1)
    for tt = 1:size(CSIT{ms},2)
        for hh = 1:size(x_corr.parameters.filters,1)

            temp = eegfilt(scale.*(data.lfp{7,1}(ii,:,CSIT{ms}(1,tt))),parameters.decimated_srate,x_corr.parameters.filters(hh,1),x_corr.parameters.filters(hh,2)); % just filtering
            x_corr.amplitude{2,hh}(ii,:,tt) = abs(hilbert(temp))-mean(abs(hilbert(temp)),2); % getting the amplitude envelope and subtracted the mean to eliminate each signal's DC component
%           x_corr.amplitude{2,hh}(ii,:) = temp;

        end
    end

end



% ITI
for ii = 1:size(data.lfp{8,1},1)
    for tt = 1:size(CSIT{ms},2)
        for hh = 1:size(x_corr.parameters.filters,1)

            temp = eegfilt(scale.*(data.lfp{8,1}(ii,:,CSIT{ms}(1,tt))),parameters.decimated_srate,x_corr.parameters.filters(hh,1),x_corr.parameters.filters(hh,2)); % just filtering
            x_corr.amplitude{3,hh}(ii,:,tt) = abs(hilbert(temp))-mean(abs(hilbert(temp)),2); % getting the amplitude envelope and subtracted the mean to eliminate each signal's DC component
%            x_corr.amplitude{3,hh}(ii,:) = temp;

        end
    end

end


clear('ii','hh','tt','temp','scale')



%% Compute Cross-correlation

% Cell colums -> amplitude filtered signals
% Cell Rows
%    . 1 -> baseline
%    . 2 -> CS-Tone
%    . 3 -> ITI

% In each cell
% First  dimention --> Combinations
%  .Row 1 -> mPFC PL <--> mPFC IL
%  .Row 2 -> mPFC PL <--> dHPC
%  .Row 3 -> mPFC IL <--> dHPC

% Second dimention --> Lags
% Third  dimention --> Time


% all possible combinations between all channels:
x_corr.parameters.combinations  = nchoosek(1:size(x_corr.amplitude{1,1},1),2);

% Compute how many time points are in one cycle, and limit corr/xcov to this lag
% Number of lags need to be equal or less than one cycle

% Parameters
window_length1 = 1.9980 %6*(1/3);
window_length2 = 1.0020 %*6*(1/6);
window_length3 = 0.8571 %6*(1/7);

x_corr.parameters.window_length = [window_length1 window_length2 window_length3];
x_corr.parameters.nlags   =  250;
x_corr.parameters.overlap =  50;
x_corr.parameters.method  = 'coeff';

fp = 0;

% Initialze
x_corr.result_full                    = cell(size(x_corr.amplitude));
x_corr.result_full_lag_peak           = cell(size(x_corr.amplitude));
x_corr.result_full_lag_peak_idx       = cell(size(x_corr.amplitude));
x_corr.result                         = cell(size(x_corr.amplitude));
x_corr.result_lag_peak                = cell(size(x_corr.amplitude));



for hh = 1:size(x_corr.amplitude,1)
    for ii = 1:size(x_corr.amplitude,2)
        for jj = 1:size(x_corr.parameters.combinations,1)

            if hh == 1 % condition for baseline

                signal1 = x_corr.amplitude{hh,ii}(x_corr.parameters.combinations(jj,1),:);
                signal2 = x_corr.amplitude{hh,ii}(x_corr.parameters.combinations(jj,2),:);

                [x_corr.result_full{hh,ii}(jj,:), x_corr.result_full_lag_peak{hh,ii}(jj,:), x_corr.result_full_lag_peak_idx{hh,ii}(jj,:), x_corr.result{hh,ii}(jj,:,:),x_corr.result_lag_peak{hh,ii}(jj,:),x_corr.lags] = ...
                    x_corr_overtime(signal1, signal2, x_corr.parameters.nlags, parameters.decimated_srate, x_corr.parameters.window_length.*10, x_corr.parameters.overlap, x_corr.parameters.method,fp); % coeff return correlation coeficients


            else % condition for CS and ITI

                for tt = 1:size(x_corr.amplitude{hh,ii},3)

                    signal1 = x_corr.amplitude{hh,ii}(x_corr.parameters.combinations(jj,1),:,tt);
                    signal2 = x_corr.amplitude{hh,ii}(x_corr.parameters.combinations(jj,2),:,tt);


                    [x_corr.result_full{hh,ii}(jj,:,tt), x_corr.result_full_lag_peak{hh,ii}(jj,tt), x_corr.result_full_lag_peak_idx{hh,ii}(jj,tt), x_corr.result{hh,ii}(jj,:,:,tt), x_corr.result_lag_peak{hh,ii}(jj,:,tt)] = ...
                        x_corr_overtime(signal1, signal2, x_corr.parameters.nlags, parameters.decimated_srate, x_corr.parameters.window_length(1,ii), x_corr.parameters.overlap, x_corr.parameters.method,fp); % coeff return correlation coeficients

                end
            end
        end
    end
end


clear('signal1','signal2','window_length1','window_length2','window_length3','fp','hh','ii','jj','tt')

%% Normalization and trials average

x_corr.result_full_trials      = [];
x_corr.result_full_trials_peak = [];


for jj = 1:size(x_corr.result_full,1)
    for ii = 1:size(x_corr.result_full,2)

        % Averaging trials

        if jj == 1 % baseline condition

            x_corr.result_full_trials{jj,ii}                   = x_corr.result_full{jj,ii};              % baseline cell will in 1,1 of course....but just to keep the organization
            x_corr.result_full_trials_peak{jj,ii}              = x_corr.result_full_lag_peak{jj,ii};
            x_corr.result_full_trials_mean{jj,ii}              = mean(x_corr.result_full_trials{jj,ii},3); 

            [x_corr.result_full_trials_r_mean{jj,ii},x_corr.result_full_trials_peak_mean_idx{jj,ii}] = max(x_corr.result_full_trials_mean{jj,ii},[],2);

            x_corr.result_full_trials_peak_mean{jj,ii}         = x_corr.lags(:,x_corr.result_full_trials_peak_mean_idx{jj,ii});

         

        else

            x_corr.result_full_trials{jj,ii}                   = x_corr.result_full{jj,ii};                           
            x_corr.result_full_trials_peak{jj,ii}              = x_corr.result_full_lag_peak{jj,ii};
            x_corr.result_full_trials_mean{jj,ii}              = mean(x_corr.result_full_trials{jj,ii},3); 

            [x_corr.result_full_trials_r_mean{jj,ii},x_corr.result_full_trials_peak_mean_idx{jj,ii}] = max(x_corr.result_full_trials_mean{jj,ii},[],2);
            
            x_corr.result_full_trials_peak_mean{jj,ii}         = x_corr.lags(:,x_corr.result_full_trials_peak_mean_idx{jj,ii});
                        
        end

    end
end


% Peaks over time

x_corr.result_full_trials_overtime{jj,ii}    = [];
x_corr.result_full_trials_peak_overtime      = [];
x_corr.result_full_trials_peak_overtime_idx  = [];

for jj = 1:size(x_corr.result_full,1)
    for ii = 1:size(x_corr.result_full,2)

        if jj == 1 % baseline condition
            x_corr.result_full_trials_overtime{jj,ii}              =  normalize(x_corr.result{jj,ii},3,'range');
            [~,x_corr.result_full_trials_peak_overtime_idx{jj,ii}] =  max(x_corr.result_full_trials_overtime{jj,ii},[],2);
            x_corr.result_full_trials_peak_overtime{jj,ii}         =  x_corr.lags(x_corr.result_full_trials_peak_overtime_idx{jj,ii});
            x_corr.result_full_trials_peak_overtime_mean{jj,ii}     =  mean(x_corr.result_full_trials_peak_overtime{jj,ii},2);

        else
            x_corr.result_full_trials_overtime{jj,ii}              =  normalize(mean(x_corr.result{jj,ii},4),3,'range');
            [~,x_corr.result_full_trials_peak_overtime_idx{jj,ii}] =  max( x_corr.result_full_trials_overtime{jj,ii},[],3);
            x_corr.result_full_trials_peak_overtime{jj,ii}         =  x_corr.lags(x_corr.result_full_trials_peak_overtime_idx{jj,ii});
            x_corr.result_full_trials_peak_overtime_mean{jj,ii}     =  mean(x_corr.result_full_trials_peak_overtime{jj,ii},2);


        end
    end
end


%% Plot the results

% choose filter
ff = 2;

%xlim([-x_corr.parameters.nlags x_corr.parameters.nlags])

for jj = 1:3

    figure(jj);
    set(gcf,'color','w');
    sc = [1,1,800,1200];
    set(gcf, 'Position', sc);

    title_ = [{'mPFC PL lead <-----> mPFC IL lead'};{'mPFC PL lead <-----> dHPC lead'};{'mPFC IL lead <-----> dHPC lead'}];
    sgtitle([{'Cross-Correlation Over Time '};title_{jj}],'FontSize',14);


    % Baseline
    subplot(5,2,1);
    plot(x_corr.lags,x_corr.result_full{1,ff}(jj,:),'LineWidth',2,'color',[.6 .6 .6]);
    hold on
    plot(x_corr.lags(x_corr.result_full_lag_peak_idx{1,ff}(jj,1)),x_corr.result_full{1,ff}(jj,x_corr.result_full_lag_peak_idx{1,ff}(jj,1)),'ro','MarkerEdgeColor','r','MarkerFaceColor','r');
    title({['lag = ', num2str(x_corr.result_full_lag_peak{1,ff}(jj,1)),' ms'];[]})
    ylabel([{'Normalized'};{'Crosscorrelation'}]);
    xlim([-x_corr.parameters.nlags x_corr.parameters.nlags])
    ylim([min(x_corr.result_full{1,ff}(jj,:)) max(x_corr.result_full{1,ff}(jj,:))])
    xline(0,'k--','LineWidth', 2)
    box off
    % legend('','peak lag','FontSize',10,'location', 'east')
    % legend('boxoff')

    subplot(5,2,3)
    histogram(x_corr.result_lag_peak{1,ff}(jj,:),'BinWidth',30,'FaceColor',[.6 .6 .6]);
    ylabel('Count');
    % a = fitdist(xcorr_result_lag_peak','normal');
    % hold on
    % plot(a)
    xlim([-x_corr.parameters.nlags x_corr.parameters.nlags])
    box off
    %ylim([0 10])
    xline(0,'k--','LineWidth', 2)

    subplot(5,2,[5 7 9])
    imagesc(x_corr.lags, (1:size(x_corr.result{1,ff},2)), normalize(squeeze(x_corr.result{1,ff}(jj,:,:)),2,'range'));
    %contourf(x_corr.lags,(1:size(x_corr.result{1,ff},2)),squeeze(x_corr.result{1,ff}(jj,:,:)),80,'linecolor','none');
    xlim([-x_corr.parameters.nlags x_corr.parameters.nlags])
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
    r=plot(x_corr.result_lag_peak{1,ff}(jj,:),(1:size(x_corr.result{1,ff},2)),'wo','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',4);
    xlim([-x_corr.parameters.nlags x_corr.parameters.nlags])
    ylim([0.5 size(x_corr.result{1,ff},2)+.5])
    a = gca; % Get axis
    a.YColor = 'w';
    a.YTick = [];
    % legend(r,'peak','location','southoutside','FontSize',14)
    % legend('boxoff')


    % CS-Tone. First 10 trials

    subplot(5,2,2)
    plot(x_corr.lags,x_corr.result_full_trials_mean{2,ff}(jj,:),'LineWidth',2,'color',[.6 .6 .6]);
    hold on
    plot(x_corr.lags(x_corr.result_full_trials_peak_mean_idx{2,ff}(jj,1)),x_corr.result_full_trials_mean{2,ff}(jj,x_corr.result_full_trials_peak_mean_idx{2,ff}(jj,1)),'ro','MarkerEdgeColor','r','MarkerFaceColor','r');
    title({['lag = ', num2str(x_corr.result_full_trials_peak_mean{2,ff}(1,jj)),' ms'];[]})
    ylabel([{'Normalized'};{'Crosscorrelation'}]);
    xlim([-x_corr.parameters.nlags x_corr.parameters.nlags])
    ylim([min(x_corr.result_full_trials_mean {2,ff}(jj,:)) max(x_corr.result_full_trials_mean {2,ff}(jj,:))])
    xline(0,'k--','LineWidth', 2)
    box off
    % legend('','peak lag','FontSize',10,'location', 'east')
    % legend('boxoff')

    subplot(5,2,4)
    histogram(x_corr.result_full_trials_peak_overtime{2,ff}(jj,:),'BinWidth',30,'FaceColor',[.6 .6 .6]);
    ylabel('Count');
    % a = fitdist(xcorr_result_lag_peak','normal');
    % hold on
    % plot(a)
    xlim([-x_corr.parameters.nlags x_corr.parameters.nlags])
    box off
    %ylim([0 10])
    xline(0,'k--','LineWidth', 2)

    subplot(5,2,[6 8 10])
    imagesc(x_corr.lags, (1:size(x_corr.result{2,ff},2)), squeeze(x_corr.result_full_trials_overtime{2,ff}(jj,:,:)));
    %contourf(x_corr.lags,(1:size(x_corr.result{2,ff},2)), mean(squeeze(x_corr.result{2,ff}(jj,:,:,1:size(CSIT{ms},2)/2)),3),80,'linecolor','none');

    xlim([-x_corr.parameters.nlags x_corr.parameters.nlags])
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
    r=plot(x_corr.result_full_trials_peak_overtime{2,ff}(jj,:),(1:size(x_corr.result{2,ff},2)),'wo','MarkerEdgeColor','k','MarkerFaceColor','w');
    xlim([-x_corr.parameters.nlags x_corr.parameters.nlags])
    ylim([0.5 size(x_corr.result{2,ff},2)+.5])
    a = gca; % Get axis
    a.YColor = 'w';
    a.YTick = [];
    % legend(r,'peak','location','southoutside','FontSize',14)
    % legend('boxoff')

end

% linkaxes([b c d],'xy');
%b(1).YLim  = [0 15];
% b(1).XLim  = [-250 250];

%% Save Figures

newStr = id(1:end-8);

path = '/Users/flavio/Desktop';
%path = files.FilesLoaded{1, 1}.folder;

name = strcat(path,'/',newStr,'_x_corr_retrieval');

%get all figures
figlist = findobj(allchild(0), 'flat', 'Type', 'figure');

%name_figs = {'_3_5Hz','_6_8Hz','_7_10Hz'};
name_figs = {'_PL_IL','_PL_dHPC','_IL_dHPC'};

set(gcf,'renderer','Painters')

% Loop through figure
for ii = 1:numel(figlist)
    name_loop = strcat(name,name_figs{ii});
    FigHandle = figlist(ii);
    saveas(FigHandle,name_loop,'png')
    exportgraphics(gcf,strcat(name_loop,'.eps'),'Resolution', 300)

end

close all
clear('name','newStr','path')

%% Save data

%newStr = regexprep(files.id.name,'.mat','_');
%newStr = files.id(ms).name(1:end-8);
newStr = id(1:end-8);

path = '/Users/flavio/Desktop';
%path = files.FilesLoaded{1, 1}.folder;

name = strcat(path,'/',newStr,'_x_corr_retrieval.mat');

% save data
save(name,'x_corr','-v7.3')

clear('name','newStr','path')

%% last update 19/06/2024
%  listening:
