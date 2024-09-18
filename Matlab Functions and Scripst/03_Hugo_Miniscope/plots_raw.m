

data.stats.raw{1,1} = data.raw{1,1};


idx = 2:2:20

% consolidation
for ii = 1:length(idx)
    data.stats.raw{2,ii}(:,:,jj) = cell2mat(temp(1,jj));
end

% Extinction 10 first ITI

idx = 4:2:22;

for ii = 1:length(idx)
    data.stats.raw{3,1}(:,:,ii) = data.raw{3,idx(ii)};
end

% Extinction 10 last ITI
idx = 74:2:92;

for ii = 1:length(idx)
    data.stats.raw{3,2}(:,:,ii) = data.raw{3,idx(ii)};
end


% Retrieval 10 first ITI
idx = 4:2:22;

for ii = 1:length(idx)
    data.stats.raw{4,1}(:,:,ii) = data.raw{4,idx(ii)}(1:145,:);
end

% Retrieval 10 last ITI
idx = 74:2:92;

for ii = 1:length(idx)
    data.stats.raw{4,2}(:,:,ii) = data.raw{4,idx(ii)}(1:145,:);
end


%% Plot

data.stats.raw_mean =[];

temp = cell2mat(data.stats.raw(1,1));
data.stats.raw_mean{1,1} = mean(temp(1:295,:),2,'omitnan');

temp = cell2mat(data.stats.raw(2,2));
data.stats.raw_mean{2,1} = mean(mean(temp,3,'omitnan'),2,'omitnan');

temp = cell2mat(data.stats.raw(3,1));
data.stats.raw_mean{3,1} = squeeze(mean(temp,2,'omitnan'));

temp = cell2mat(data.stats.raw(3,2));
data.stats.raw_mean{3,2} = squeeze(mean(temp,2,'omitnan'));

temp = cell2mat(data.stats.raw(4,1));
data.stats.raw_mean{4,1} = squeeze(mean(temp,2,'omitnan'));

temp = cell2mat(data.stats.raw(4,2));
data.stats.raw_mean{4,2} = squeeze(mean(temp,2,'omitnan'));


timev_baseline      = linspace(1,300,length(data.stats.raw_mean{1,1}));
timev_concolidation = linspace(1,300,length(data.stats.raw_mean{2,1}));
timev_extinction    = linspace(1,30,length(data.stats.raw_mean{3,1}));
timev_retrieval     = linspace(1,30,length(data.stats.raw_mean{4,1}));


%%
figure
subplot(5,10,[1 5])
plot(timev_baseline,data.stats.raw_mean{1,1})
ylim([0 10])

subplot(5,10,[6 10])
plot(timev_concolidation,data.stats.raw_mean{2,1})
ylim([0 10])

for ii = 1:size(data.stats.raw_mean{3,1},2)
    subplot(5,10,10+ii)
    plot(timev_extinction,data.stats.raw_mean{3,1}(:,ii))
    ylim([0 10])
end

for ii = 1:size(data.stats.raw_mean{3,2},2)
    subplot(5,10,20+ii)
    plot(timev_retrieval,data.stats.raw_mean{3,2}(:,ii),'r')
    ylim([0 10])
end

for ii = 1:size(data.stats.raw_mean{4,1},2)
    subplot(5,10,30+ii)
    plot(timev_extinction,data.stats.raw_mean{4,1}(:,ii))
    ylim([0 10])
end

for ii = 1:size(data.stats.raw_mean{4,2},2)
    subplot(5,10,40+ii)
    plot(timev_retrieval,data.stats.raw_mean{4,2}(:,ii),'r')
    ylim([0 10])
end

%% Acitivity
data.stats.raw{1,1} = data.raw{1,1};


idx = 2:2:20

% consolidation
for ii = 1:length(idx)
    data.stats.raw{2,ii}(:,:,jj) = cell2mat(temp(1,jj));
end

% Extinction 10 first ITI

idx = 4:2:22;

for ii = 1:length(idx)
    data.stats.raw{3,1}(:,:,ii) = data.raw{3,idx(ii)};
end

% Extinction 10 last ITI
idx = 74:2:92;

for ii = 1:length(idx)
    data.stats.raw{3,2}(:,:,ii) = data.raw{3,idx(ii)};
end


% Retrieval 10 first ITI
idx = 4:2:22;

for ii = 1:length(idx)
    data.stats.raw{4,1}(:,:,ii) = data.raw{4,idx(ii)}(1:145,:);
end

% Retrieval 10 last ITI
idx = 74:2:92;

for ii = 1:length(idx)
    data.stats.raw{4,2}(:,:,ii) = data.raw{4,idx(ii)}(1:145,:);
end


%% Plot

data.stats.raw_mean =[];

temp = cell2mat(data.stats.raw(1,1));
data.stats.raw_mean{1,1} = mean(temp(1:295,:),2,'omitnan');

temp = cell2mat(data.stats.raw(2,2));
data.stats.raw_mean{2,1} = mean(mean(temp,3,'omitnan'),2,'omitnan');

temp = cell2mat(data.stats.raw(3,1));
data.stats.raw_mean{3,1} = squeeze(mean(temp,2,'omitnan'));

temp = cell2mat(data.stats.raw(3,2));
data.stats.raw_mean{3,2} = squeeze(mean(temp,2,'omitnan'));

temp = cell2mat(data.stats.raw(4,1));
data.stats.raw_mean{4,1} = squeeze(mean(temp,2,'omitnan'));

temp = cell2mat(data.stats.raw(4,2));
data.stats.raw_mean{4,2} = squeeze(mean(temp,2,'omitnan'));


timev_baseline      = linspace(1,300,length(data.stats.raw_mean{1,1}));
timev_concolidation = linspace(1,300,length(data.stats.raw_mean{2,1}));
timev_extinction    = linspace(1,30,length(data.stats.raw_mean{3,1}));
timev_retrieval     = linspace(1,30,length(data.stats.raw_mean{4,1}));


%%
figure
subplot(5,10,[1 5])
plot(timev_baseline,data.stats.raw_mean{1,1})
ylim([0 10])

subplot(5,10,[6 10])
plot(timev_concolidation,data.stats.raw_mean{2,1})
ylim([0 10])

for ii = 1:size(data.stats.raw_mean{3,1},2)
    subplot(5,10,10+ii)
    plot(timev_extinction,data.stats.raw_mean{3,1}(:,ii))
    ylim([0 10])
end

for ii = 1:size(data.stats.raw_mean{3,2},2)
    subplot(5,10,20+ii)
    plot(timev_retrieval,data.stats.raw_mean{3,2}(:,ii),'r')
    ylim([0 10])
end

for ii = 1:size(data.stats.raw_mean{4,1},2)
    subplot(5,10,30+ii)
    plot(timev_extinction,data.stats.raw_mean{4,1}(:,ii))
    ylim([0 10])
end

for ii = 1:size(data.stats.raw_mean{4,2},2)
    subplot(5,10,40+ii)
    plot(timev_retrieval,data.stats.raw_mean{4,2}(:,ii),'r')
    ylim([0 10])
end




