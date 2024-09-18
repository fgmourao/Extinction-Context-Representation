%% Analysis of cellular activity. Calcium imaging experiment. Miniscope.
%  Pre-Processing - Sorting Session, events and cells

%  After preprocessing, the data.activity from both sessions_activity were merged and saved in the same spreadsheet.
%  the table includes time (seconds) where the event occurs, the cell ID that produced the event, and the amplitude (value) of each event.

%  The MedPC trigger starts the recording, so 0s on the miniscope recording is also 0s on the freezing output


% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  02/2024
% Last update: 03/2024

%% Stats

data.stats.all_cells{1,1} = data.activity_sorted_baseline(5,:);         % baseline
data.stats.all_cells{2,1} = data.activity_sorted_consolidation(5,:);    % consolidation
data.stats.all_cells{3,1} = data.activity_sorted_extinction(5,:);       % extinction
data.stats.all_cells{4,1} = data.activity_sorted_retrieval(5,:);        % retreival

for ii = 1:size(data.stats.all_cells,1)
    for jj = 1:size(data.stats.all_cells{ii,1},2)
        data.stats.all_cells{ii,2}{1,jj} = sum(data.stats.all_cells{ii,1}{1,jj},2); % Total fire for each cell. That is, sum over time
        data.stats.all_cells{ii,3}{1,jj} = sum(data.stats.all_cells{ii,1}{1,jj},1); % Total fire over time. That is, each cell was summed


    end
end


% Histograms and bin edges
% Baseline
for jj = 1:size(data.stats.all_cells{1,1},2)
    
    if jj == 1
        bins = 300; % 5-s time epoch for entire session
    else
        bins = 30; % 2-s time epoch for each block
    end

    data.stats.all_cells{1,4}{1,jj} = linspace(0,round(size(data.stats.all_cells{1, 1}{1, jj},2)./parameters.sampling_rate),bins); % Edges
    data.stats.all_cells{1,5}{1,jj} = zeros(1,length(data.stats.all_cells{1,4}{1,jj})); % Counts

    % Histogram counts over time
    for mm = 1:size(data.stats.all_cells{1, 1}{1, jj},1)
        timev_raw = data.raw{1,jj}(:,1)';
        data.stats.all_cells{1,5}{1,jj}= data.stats.all_cells{1,5}{1,jj} + histc(timev_raw(data.stats.all_cells{1, 1}{1, jj}(mm,:)==1),data.stats.all_cells{1,4}{1,jj});
    end

    % Fit Histogram count
    data.stats.all_cells{1,6}{1,jj} = smooth(data.stats.all_cells{1,4}{1,jj} ,data.stats.all_cells{1,5}{1,jj},0.4,'rloess');

    % Fit probability distribution over Histogram counts
    data.stats.all_cells{1,7}{1,jj} = fitdist(data.stats.all_cells{1,5}{1,jj}(data.stats.all_cells{1,5}{1,jj}>0)','Lognormal');

    % Histogram counts to probability distribution
    [data.stats.all_cells{1,9}{1,jj},data.stats.all_cells{1,8}{1,jj}] = histcounts(data.stats.all_cells{1,5}{1,jj});

    % Probability density function
    data.stats.all_cells{1,10}{1,jj} = pdf(data.stats.all_cells{1,7}{1,jj},data.stats.all_cells{1,8}{1,jj});
    
    % Cumulative density function
    data.stats.all_cells{1,11}{1,jj} = cdf(data.stats.all_cells{1,7}{1,jj},data.stats.all_cells{1,8}{1,jj});

end

% Histograms and bin edges
% Consolidation
for jj = 1:size(data.stats.all_cells{2,1},2)
    
    if jj == 1
        bins = 1200; % 5-s time epoch for entire session
    else
        bins = 300; % 2-s time epoch for each block
    end

    data.stats.all_cells{2,4}{1,jj} = linspace(0,round(size(data.stats.all_cells{2, 1}{1, jj},2)./parameters.sampling_rate),bins); % Edges
    data.stats.all_cells{2,5}{1,jj} = zeros(1,length(data.stats.all_cells{2,4}{1,jj})); % Counts

    % Histogram counts over time
    for mm = 1:size(data.stats.all_cells{2, 1}{1, jj},1)
        timev_raw = data.raw{2,jj}(:,1)';
        data.stats.all_cells{2,5}{1,jj}= data.stats.all_cells{2,5}{1,jj} + histc(timev_raw(data.stats.all_cells{2, 1}{1, jj}(mm,:)==1),data.stats.all_cells{2,4}{1,jj});
    end

    % Fit Histogram count
    data.stats.all_cells{2,6}{1,jj} = smooth(data.stats.all_cells{2,4}{1,jj} ,data.stats.all_cells{2,5}{1,jj},0.4,'rloess');

    % Fit probability distribution over Histogram counts
    data.stats.all_cells{2,7}{1,jj} = fitdist(data.stats.all_cells{2,5}{1,jj}(data.stats.all_cells{2,5}{1,jj}>0)','Lognormal');

    % Histogram counts to probability distribution
    [data.stats.all_cells{2,9}{1,jj},data.stats.all_cells{2,8}{1,jj}] = histcounts(data.stats.all_cells{2,5}{1,jj});

    % Probability density function
    data.stats.all_cells{2,10}{1,jj} = pdf(data.stats.all_cells{2,7}{1,jj},data.stats.all_cells{2,8}{1,jj});
    
    % Cumulative density function
    data.stats.all_cells{2,11}{1,jj} = cdf(data.stats.all_cells{2,7}{1,jj},data.stats.all_cells{2,8}{1,jj});

end


% Histograms and bin edges
% Extinction
for jj = 1:size(data.stats.all_cells{3,1},2)
    
    if jj == 1
        bins = 60; % 5-s time epoch for entire session
    else
        bins = 30; % 2-s time epoch for each block
    end

    data.stats.all_cells{3,4}{1,jj} = linspace(0,round(size(data.stats.all_cells{3, 1}{1, jj},2)./parameters.sampling_rate),bins); % Edges
    data.stats.all_cells{3,5}{1,jj} = zeros(1,length(data.stats.all_cells{3,4}{1,jj})); % Counts

    % Histogram counts over time
    for mm = 1:size(data.stats.all_cells{3, 1}{1, jj},1)
        timev_raw = data.raw{3,jj}(:,1)';
        data.stats.all_cells{3,5}{1,jj}= data.stats.all_cells{3,5}{1,jj} + histc(timev_raw(data.stats.all_cells{3, 1}{1, jj}(mm,:)==1),data.stats.all_cells{3,4}{1,jj});
    end

%     % Fit Histogram count
    data.stats.all_cells{3,6}{1,jj} = smooth(data.stats.all_cells{3,4}{1,jj} ,data.stats.all_cells{3,5}{1,jj},0.4,'rloess');
% 
%     % Fit probability distribution over Histogram counts
%     data.stats.all_cells{3,7}{1,jj} = fitdist(data.stats.all_cells{3,5}{1,jj}(data.stats.all_cells{3,5}{1,jj}>0)','Lognormal');
% 
%     % Histogram counts to probability distribution
%     [data.stats.all_cells{3,9}{1,jj},data.stats.all_cells{3,8}{1,jj}] = histcounts(data.stats.all_cells{3,5}{1,jj});
% 
%     % Probability density function
%     data.stats.all_cells{3,10}{1,jj} = pdf(data.stats.all_cells{3,7}{1,jj},data.stats.all_cells{3,8}{1,jj});
%     
%     % Cumulative density function
%     data.stats.all_cells{3,11}{1,jj} = cdf(data.stats.all_cells{3,7}{1,jj},data.stats.all_cells{3,8}{1,jj});

end


% Histograms and bin edges
% Retrival
for jj = 1:size(data.stats.all_cells{4,1},2)
    
    if jj == 1
        bins = 60; % 5-s time epoch for entire session
    else
        bins = 30; % 2-s time epoch for each block
    end

    data.stats.all_cells{4,4}{1,jj} = linspace(0,round(size(data.stats.all_cells{4, 1}{1, jj},2)./parameters.sampling_rate),bins); % Edges
    data.stats.all_cells{4,5}{1,jj} = zeros(1,length(data.stats.all_cells{4,4}{1,jj})); % Counts

    % Histogram counts over time
    for mm = 1:size(data.stats.all_cells{4, 1}{1, jj},1)
        timev_raw = data.raw{4,jj}(:,1)';
        data.stats.all_cells{4,5}{1,jj}= data.stats.all_cells{4,5}{1,jj} + histc(timev_raw(data.stats.all_cells{4, 1}{1, jj}(mm,:)==1),data.stats.all_cells{4,4}{1,jj});
    end

    % Fit Histogram count
    data.stats.all_cells{4,6}{1,jj} = smooth(data.stats.all_cells{4,4}{1,jj} ,data.stats.all_cells{4,5}{1,jj},0.4,'rloess');

%     % Fit probability distribution over Histogram counts
%     data.stats.all_cells{4,7}{1,jj} = fitdist(data.stats.all_cells{4,5}{1,jj}(data.stats.all_cells{4,5}{1,jj}>0)','Lognormal');
% 
%     % Histogram counts to probability distribution
%     [data.stats.all_cells{4,9}{1,jj},data.stats.all_cells{4,8}{1,jj}] = histcounts(data.stats.all_cells{4,5}{1,jj});
% 
%     % Probability density function
%     data.stats.all_cells{4,10}{1,jj} = pdf(data.stats.all_cells{4,7}{1,jj},data.stats.all_cells{4,8}{1,jj});
%     
%     % Cumulative density function
%     data.stats.all_cells{4,11}{1,jj} = cdf(data.stats.all_cells{4,7}{1,jj},data.stats.all_cells{4,8}{1,jj});

end

clear ('bins','ii','timev_raw')

%% Plots

% choose data
ds = 3; % session
db = 1; % block

% Choose Trial
% data 1 - raw signal intensity
data2plot_raw = data.raw{ds, db}(:,2:end)';
rowsToKepp = all(~isnan(data2plot_raw), 2);
data2plot_raw = data2plot_raw(rowsToKepp,:);
%data2plot_raw(any(isnan(data2plot_raw), 2), :) = [];

[data2plot_raw_normalized, ~]=mapminmax(data2plot_raw,0,1);
timev_raw = linspace(0,size(data2plot_raw,2)/parameters.sampling_rate,size(data2plot_raw,2));

% data 2 - raster plots
data2plot_raster = data.stats.all_cells{ds, 1}{1,db}(rowsToKepp,:);
data2plot_raster(data2plot_raster == 0) = nan;
timev_raster = linspace(0,size(data2plot_raster,2)/parameters.sampling_rate,size(data2plot_raster,2));


% data 3 = histogram activity over time
data2plot_hist_edges  = data.stats.all_cells{ds,4}{1,db};
data2plot_hist_counts = data.stats.all_cells{ds,5}{1,db};
data2plot_hist_fit = data.stats.all_cells{ds,6}{1,db};

% data 3 - probability
data2plot_prob_pd = data.stats.all_cells{ds,7}{1,db};
data2plot_prob_pedges = data.stats.all_cells{ds,8}{1,db};
data2plot_prob_pcounts = data.stats.all_cells{ds,9}{1,db};
data2plot_prob_pdf = data.stats.all_cells{ds,10}{1,db};
data2plot_prob_cdf = data.stats.all_cells{ds,11}{1,db};



%Cells
cells_ = 1:size(data2plot_raw_normalized,1);

% factor
factor_1 = (cells_)'*1;
factor_2 = (cells_)'*.5;

% Title
t = {'Baseline','Consolidation','Extinction','retrieval'};

figure
set(gcf,'color','w');
sc = [1,1,960,1200];
set(gcf, 'Position', sc);

% Raw signal waveforms
subplot(4,4,[1 2 5 6])
sgtitle(t(ds))
hold on

r = plot(timev_raw, bsxfun(@plus, data2plot_raw_normalized(cells_,:), factor_1),'Color','[0.3, 0.3, 0.3]','linew',1);
a = gca;
a.YColor = 'w';
a.YTick = [];
%ylabel('\fontsize{12}Cells (1 -> 57)');
a.XLim = [0  timev_raw(end)];
box off


% Raster
subplot(4,4,[3 4 7 8])

r = plot(timev_raster, bsxfun(@plus, data2plot_raster(cells_,:), factor_2),"|",'Color','[0.3, 0.3, 0.3]','linew',1,'MarkerSize',5);
a = gca;
a.YColor = 'w';
a.YTick = [];
a.XLim = [0  timev_raw(end)];
box off


% Raw signal intensity
subplot(4,4,[9 10 13 14])

contourf(timev_raw,cells_,data2plot_raw(cells_,:),80,'linecolor','none');
a = gca;
%a.YColor = 'w';
a.YTick = [];
%xlim([150 165])
clim([0 max(max(data2plot_raw(cells_,:)))-20])
xlabel('\fontsize{12}Time (s)');
%ylabel('\fontsize{12}Cells (1 -> 57)');

c = colorbar;
set(get(c,'ylabel'),'string','\fontsize{12}\DeltaF/F','Rotation',270); 
set(c,'XTickLabel',{'Low',' ',' ',' ',' ',' ',' ',' ',' ',' ','High'}); 
%view(0,90)
py_path = "~/anaconda3/envs/Python_3_10/bin/python";
Py_map = getPyPlot_cMap('viridis', [], [], py_path);
colormap(Py_map)


% PEHS
subplot(4,4,[11 12])

b1 = bar(data2plot_hist_edges,data2plot_hist_counts);
b1(1).FaceColor = [.6 .6 .6,];
b1(1).LineWidth = .5;
%xlim([0 300])
%ylim([0 100])
xlabel('\fontsize{12}Time (s)');
ylabel('\fontsize{12}Frequency');
hold on
plot(data2plot_hist_edges,data2plot_hist_fit,'LineWidth',2,'color', [0.8500 0.3250 0.0980]);
legend('', 'rloess fit')
legend('boxoff')

% PDF
subplot(4,4,15)

pp = plot(data2plot_prob_pd)
pp(1).LineWidth = 2;
pp(2).FaceColor = 'w';
pp(2).LineWidth = 1;
pp(2).EdgeColor ='k';
xlabel('Activity per bin','FontSize',12), ylabel(' Probability','FontSize',12)
xlim([pp(2).BinEdges(1)-1 pp(2).BinEdges(end)+1 ])
legend('', 'lognormal')
legend('boxoff')

% CDF
subplot(4,4,16)
plot(data2plot_prob_pedges,data2plot_prob_cdf,'LineWidth',2,'color',[0.8500 0.3250 0.0980] )
xlabel('Activity per bin','FontSize',12), ylabel('Cumulative Probability','FontSize',12)


%% 
% Flavio Mourao. Last update 23/04/18. 16.10am. David Wilson Library. Leicester Uk
% listening: Mogwai - May nothing but hapiness come through your door 