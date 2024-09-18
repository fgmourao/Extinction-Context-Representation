
%% Correlation and Covariance Matrices betwwen channels

% - Performs Correlation and Covariance Matrices betwwen channels 
%   considering the trial periods
% - Power Spearman's correlation betwwen channels
%   considering the trial periods

% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  12/2023
% Last update:

%%
% First extract the data with: Extracting_LFPs_and_events.m
% ... then organize the data with the script: Pre_processing.m



[Zica, W, T, mu] = fastICA(data.lfp{5, 1},32,'kurtosis',1);


%% check all ICA components

figure

for ii = 1:size(Zica,1)

    subplot(4,8,ii)
    plot(Zica(ii,:));
    title(num2str(ii))

end


%% check all channels

figure

for ii = 1:size(data.lfp{5, 1},1)

    subplot(4,8,ii)
    plot(data.lfp{5, 1}(ii,:));
    title(num2str(ii))

end


%% Use first channel from each probe as Reference

for ii = 5:size(data.lfp,1)

    temp1 = data.lfp{ii, 1}(1:8,:,:) - data.lfp{ii, 1}(1,:,:);    % mPFC IL
    temp2 = data.lfp{ii, 1}(9:16,:,:) - data.lfp{ii, 1}(9,:,:);   % mPFC PL
    temp3 = data.lfp{ii, 1}(17:32,:,:) - data.lfp{ii, 1}(22,:,:); % dHPC

    test_plot{ii, 1} = [temp1;temp2;temp3];

end
%%
figure

for ii = 1:size(data.lfp{5, 1},1)

    subplot(4,8,ii)
    plot(test_plot{6, 1}(ii,:));
    title(num2str(ii))
    %ylim([-500 500])

end

%% Replace data

data.lfp = test_plot;

%% 

test_plot = data.lfp;

%% Remove channels
ch_ = [1:32];
ch_2_keep = [2 3 4 13 14 18 21 24 25 27 28 29 30 31 32];

ch_(ch_2_keep) = [];

for ii = 5:size(data.lfp,1)
        data.lfp{ii,1}(ch_,:,:) = [];        

end

%%

pw = [];

% Time window
pw.full_trial.baseline_timewin    = 2048; % in ms

% Convert time window to points
pw.full_trial.baseline_timewinpnts   = hamming(round(pw.full_trial.baseline_timewin/(1000/parameters.decimated_srate)));

% nFFT
pw.full_trial.nFFT = 2^15; %4096; %2^nextpow2(pw.full_trial.baseline_timewinpnts));

% Number of overlap samples
pw.full_trial.overlap = 90;
pw.full_trial.baseline_noverlap = floor(pw.full_trial.overlap*0.01 * pw.full_trial.baseline_timewin);

% CS-Trials or freezing epochs
% Choose data from data.lfp according pre_processing.m define
% Baseline
not1 = 6;

for ii = 1:size(data.lfp{not1, 1},1)

    if ii == 1
        [pw.full_trial.Pxx{1,1}(ii,:),pw.full_trial.freq_baseline] = pwelch(test_plot{not1, 1}(ii,B_clean{ms}(1):B_clean{ms}(2)),pw.full_trial.baseline_timewinpnts,pw.full_trial.baseline_noverlap,pw.full_trial.nFFT,parameters.decimated_srate);

    else
        pw.full_trial.Pxx{1,1}(ii,:) = pwelch(test_plot{not1, 1}(ii,B_clean{ms}(1):B_clean{ms}(2)),pw.full_trial.baseline_timewinpnts,pw.full_trial.baseline_noverlap,pw.full_trial.nFFT,parameters.decimated_srate);

    end

end

not1 = 7;

for ii = 1:size(data.lfp{not1, 1},1)
    for jj = 1:size(data.lfp{not1, 1},3)

        if ii == 1
            [pw.full_trial.Pxx{2,1}(ii,:,jj),pw.full_trial.freq_baseline] = pwelch(test_plot{not1, 1}(ii,:,jj),pw.full_trial.baseline_timewinpnts,pw.full_trial.baseline_noverlap,pw.full_trial.nFFT,parameters.decimated_srate);

        else
            pw.full_trial.Pxx{2,1}(ii,:,jj) = pwelch(test_plot{not1, 1}(ii,:,jj),pw.full_trial.baseline_timewinpnts,pw.full_trial.baseline_noverlap,pw.full_trial.nFFT,parameters.decimated_srate);

        end
    end

end

not1 = 8;

for ii = 1:size(data.lfp{not1, 1},1)
    for jj = 1:size(data.lfp{not1, 1},3)

        if ii == 1
            [pw.full_trial.Pxx{3,1}(ii,:,jj),pw.full_trial.freq_baseline] = pwelch(test_plot{not1, 1}(ii,:,jj),pw.full_trial.baseline_timewinpnts,pw.full_trial.baseline_noverlap,pw.full_trial.nFFT,parameters.decimated_srate);

        else
            pw.full_trial.Pxx{3,1}(ii,:,jj) = pwelch(test_plot{not1, 1}(ii,:,jj),pw.full_trial.baseline_timewinpnts,pw.full_trial.baseline_noverlap,pw.full_trial.nFFT,parameters.decimated_srate);

        end
    end

end

%% 

theta_range = [2 12];
closestfreq_1  = find(theta_range(1)<pw.full_trial.freq_baseline & pw.full_trial.freq_baseline<theta_range(2));
freq_v = pw.full_trial.freq_baseline(closestfreq_1,1);

figure

for ii = 1:size(pw.full_trial.Pxx{1,1},1)
    subplot (4,8,ii)
    plot(freq_v, mean(pw.full_trial.Pxx{3,1}(ii,closestfreq_1,1)./sum(pw.full_trial.Pxx{3,1}(ii,closestfreq_1,1),2),3),'color',[.6, 0, 0],'linew',.5)
    xlabel('(Hz)','FontSize',9)
    %ylim([0 2000])
    %title(ch_2_keep(ii))
    title(num2str(ii))
    %ylim([0 0.01])
    xlim([2 12])

end


%% Remove channels

ch_ = [1:size(data.lfp{6, 1},1)];
ch_2_keep = [6 9 17];

ch_(ch_2_keep) = [];

for ii = 5:size(data.lfp,1)
        data.lfp{ii,1}(ch_,:,:) = [];        

end


%% remove trials

for ii = 7:size(data.lfp,1)
        data.lfp{ii,1}(:,:,11:19) = [];        

end

%% remove offset

for ii = 5:size(data.lfp,1)
    for jj = 1:size(data.lfp{ii,1},1)
        for tt = 1:size(data.lfp{ii,1},3)
            data.lfp{ii, 1}(jj,:,tt) = detrend(data.lfp{ii, 1}(jj,:,tt));
        end
    end
end

 data.lfp{9,1}(:,:,11) = [];
%% Correlation
% Choose filter band
ff = 1;

% Choose channels
ch = 1:size(data.lfp{5, 1},1);

% Choose events
ee = round(data.events{2, 1});
ee_CS = 1:10;

%ch_2_remove = [1];

figure

% baseline
data_2plot_baseline = (data.lfp{6, ff}(ch, 1:ee(1,1)-1)');
%data_2plot_baseline = (data.lfp{6, ff}(ch, 1:30000)');

%data_2plot_baseline (:,ch_2_remove) = [];

% Spearman Correlation
[correlation.cormat_base,correlation.PVAL_base] = corr(data_2plot_baseline,data_2plot_baseline,'Type','s');

subplot (2,3,1)
imagesc(correlation.cormat_base);
xlabel ('channels')
ylabel ('channels')
xticks([1:length(ch)])
yticks([1:length(ch)])

title('Baseline')
colorbar
caxis([-1 1])


% CS-Trials
for ii = 1:size(ee_CS,2)
    data_2plot_CS = data.lfp{7, ff}(ch, : ,ee_CS(ii))';

    % Spearman Correlation
    [correlation.cormat_CS,correlation.PVAL_CS] = corr( data_2plot_CS, data_2plot_CS,'Type','s');

    subplot (2,3,ii+1)
    sgtitle('Spearman Correlation')

    imagesc(correlation.cormat_CS);

    xlabel ('channels')
    ylabel ('channels')
    xticks([1:length(ch)])
    yticks([1:length(ch)])


    title(['CS-Trial: ',num2str(ee_CS(ii))])
    colorbar
    caxis([-1 1])
end

%% Save

%newStr = regexprep(files.id.name,'.mat','_');
%newStr = files.id(ms).name(1:end-8);
newStr = id(1:end-5);

path = '/Users/flavio/Desktop';
%path = files.FilesLoaded{1, 1}.folder;

name = strcat(path,'/',newStr);
% name = strcat(path,'/',newStr,'_data_retrieval');
% name = strcat(path,'/',newStr,'_data_Exposure');

save(name,'data','parameters','id','-v7.3')

clear('newStr','path')

%% Covariance

% % Choose data set
% dd = 6;
% % Choose filter band
% ff = 1;
% 
% % Choose channels
% ch = 1:29;
% 
% % Choose events
% ee = data.events{2, 1};
% 
% % Correlation
% 
% figure
% 
% % baseline
% all_channels_baseline = (data.lfp{dd, ff}(ch, 1:ee(1,1)-1)');
% 
% % Spearman Correlation
% Covariation.cormat_base = cov(all_channels_baseline);
% 
% subplot (2,3,1)
% imagesc(Covariation.cormat_base);
% 
% xlabel ('channels')
% ylabel ('channels')
% xticks([1:length(ch)])
% yticks([1:length(ch)])
% 
% title('Baseline')
% colorbar
% %caxis([0 1000])
% 
% 
% for ii = 1:size(ee,1)
%     all_channels_CS = (data.lfp{dd, ff}(ch,ee(ii,1) : ee(ii,2))');
% 
%     % Spearman Correlation
%     Covariation.cormat_CS = cov( all_channels_CS);
% 
%     subplot (2,3,ii+1)
%     sgtitle('Covariance')
% 
%     imagesc(Covariation.cormat_CS);
% 
%     xlabel ('channels')
%     ylabel ('channels')
%     xticks([1:length(ch)])
%     yticks([1:length(ch)])
% 
%     title(['CS-Trial: ',num2str(ii)])
%     colorbar
%     caxis([0 1000])
% end

% %% Exclude bad channels
% 
% % Define bad channels
% correlation.badchannels = 1;
% 
% % Exclude bad channels
% goodchannels = all_channels;
% goodchannels(:,correlation.badchannels) = [];
% 
% % Plot to check
% correlation.cormat2 = corr(goodchannels,goodchannels,'Type','s');
% 
% subplot (1,2,2)
% imagesc(correlation.cormat2);
% xlabel ('channels')
% ylabel ('channels')
% title ('Spearman Correlation only Good Channels')
% 
% clear ('ch','substrate','filter','goodchannels','all_channels')
% 
% % Exclude bad channels from data
% 
% for ii = 2:length(data.lfp(:))
% 
%     if isempty(data.lfp{ii})
%         continue
%     else
%         data.lfp{ii}(correlation.badchannels,:) = [];
%     end
% 
% end
% 
% % Correct parameters 
% parameters.nch = parameters.nch - length(correlation.badchannels);

%%

%%
temp1 = detrend(data.lfp{7, 1}(1:8,:,:) - data.lfp{7, 1}(8,:,:));    % mPFC IL
temp2 = detrend(data.lfp{7, 1}(9:16,:,:) - data.lfp{7, 1}(9,:,:));   % mPFC PL
temp3 = detrend(data.lfp{7, 1}(17:32,:,:) - data.lfp{7, 1}(17:32,:,:)); % dHPC


figure 
for ii = 1:5
    subplot(1,5,ii)
    plot(data.lfp{7, 1}(2,:,ii))
end

sg = [];
for ii = 1:5
    sg(1,:,ii) = sgolayfilt(data.lfp{7, 1}(3,:,ii),5,2001);
    subplot(1,5,ii)
    hold on
    plot(sg(1,:,ii),'LineWidth',1)
end

test = [];
for ii = 1:5
    test = data.lfp{7, 1}(3,:,ii) - sg(1,:,ii);
end


for ii = 1:size(test,3)
    [a(1,:,ii),f] = pwelch(test(1,:,ii),pw.full_trial.parameters.ITI_timewinpnts,pw.full_trial.parameters.ITI_noverlap,pw.full_trial.parameters.nFFT,parameters.decimated_srate);
end

figure
for ii = 1:5
    subplot(1,5,ii)
    plot(f,a(1,:,ii))
    xlim([2,12])
end

data.lfp{5, 1} = [];
data.lfp{5, 1} = [temp1;temp2;temp3];

%%

temp1 = detrend(data.lfp{7, 1}(1:8,:,:) - data.lfp{7, 1}(8,:,:));    % mPFC IL
temp2 = detrend(data.lfp{7, 1}(9:16,:,:) - data.lfp{7, 1}(8,:,:));   % mPFC PL
temp3 = detrend(data.lfp{7, 1}(17:32,:,:) - data.lfp{7, 1}(8,:,:)); % dHPC

data.lfp{7, 1} = [];
data.lfp{7, 1} = [temp1;temp2;temp3];

temp1 = detrend(data.lfp{9, 1}(1:8,:,:) - data.lfp{9, 1}(8,:,:));    % mPFC IL
temp2 = detrend(data.lfp{9, 1}(9:16,:,:) - data.lfp{9, 1}(8,:,:));   % mPFC PL
temp3 = detrend(data.lfp{9, 1}(17:32,:,:) - data.lfp{9, 1}(8,:,:)); % dHPC

data.lfp{9, 1} = [];
data.lfp{9, 1} = [temp1;temp2;temp3];

%% last update 18/01/2024 - 18:41
%  listening:
