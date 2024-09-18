
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

%% Channels to remove
%ch_2_remove = [1];

% Choose filter band
ff = 2;

% Choose channels
ch = 1:size(data.lfp{5, 1},1);

% Exposure

ch_2_remove = [1 4 5 8 12];

data_2plot_ = (data.lfp{5, ff}(ch, :)');
data_2plot_ (:,ch_2_remove) = [];

% Spearman Correlation
[correlation.cormat_base,correlation.PVAL_base] = corr(data_2plot_,data_2plot_,'Type','s');

imagesc(correlation.cormat_base);
xlabel ('channels')
ylabel ('channels')
xticks([1:length(ch)])
yticks([1:length(ch)])

title('Baseline')
colorbar
caxis([-1 1])

%% Extinction, Retrieval
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
%[correlation.cormat_base,correlation.PVAL_base] = corr(data_2plot_baseline,data_2plot_baseline,'Type','s');
[correlation.cormat_base] = cov(data_2plot_baseline);

%subplot (2,3,1)
imagesc(correlation.cormat_base);
xlabel ('channels')
ylabel ('channels')
xticks([1:length(ch)])
yticks([1:length(ch)])

title('Baseline')
colorbar
%caxis([-1 1])


% CS-Trials
for ii = 1:size(ee_CS,2)
    data_2plot_CS = data.lfp{7, ff}(ch, : ,ee_CS(ii))';
    %data_2plot_CS (:,ch_2_remove) = [];

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
ms = 1
%newStr = regexprep(files.id.name,'.mat','_');
newStr1 = files.id(ms).name(1:end-16);
path = files.FilesLoaded{1, 1}.folder;
%path = '/Users/flavio/Desktop';

%name = strcat('E:\Projetos 2\Flavio\Samir\Analysis\Terceiro dia\',newStr1,newStr2,'_pw_mean_Trials_allCh');
name = strcat(path,'/',newStr1,'_S_Corr_later_retrieval');

saveas(gcf,name,'png')

close all

clear('name','newStr1','path','ms')


clear('data_2plot_baseline','ii','data_2plot_CS','ee','ee_CS','ch','ff')

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

%% last update 18/01/2024 - 18:41
%  listening:
