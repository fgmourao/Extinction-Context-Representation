data.stats.activity_mean = [];

data.stats.activity_mean{1,1} = data.stats.all_cells{1, 5}{1,1};
data.stats.activity_mean{2,1} = data.stats.all_cells{2, 5}{1,1}(1,1:300);

idx = 4:2:22;
data.stats.activity_mean{3,1} = data.stats.all_cells{3, 5}(1,idx);
data.stats.activity_mean{3,1} = cat(3,data.stats.activity_mean{3, 1}{1,:});


idx = 74:2:92;
data.stats.activity_mean{3,2} = data.stats.all_cells{3, 5}(1,idx);
data.stats.activity_mean{3,2} = cat(3,data.stats.activity_mean{3, 2}{1,:});


idx = 4:2:22;
data.stats.activity_mean{4,1} = data.stats.all_cells{4, 5}(1,idx);
data.stats.activity_mean{4,1} = cat(3,data.stats.activity_mean{4, 1}{1,:});

idx = 74:2:92;
data.stats.activity_mean{4,2} = data.stats.all_cells{4, 5}(1,idx);
data.stats.activity_mean{4,2} = cat(3,data.stats.activity_mean{4, 2}{1,:});


%%

timev_baseline      = linspace(0,300,size(data.stats.activity_mean{1,1},2));
timev_concolidation = linspace(0,300,size(data.stats.activity_mean{2,1},2));
timev_extinction    = linspace(0,30,size(data.stats.activity_mean{3,1},2));
timev_retrieval     = linspace(0,30,size(data.stats.activity_mean{4,1},2));


%%
figure
set(gcf,'color','w');
sc = [1,1,960,1200];
set(gcf, 'Position', sc);


subplot(2,4,[1 2])
%plot(timev_baseline,data.stats.activity_mean{1,1})
bar(timev_baseline,data.stats.activity_mean{1,1},'FaceColor',[.6 .6 .6,])
hold on
plot(timev_baseline,smooth(data.stats.activity_mean{1,1},0.4,'rloess'),'LineWidth',2);
ylim([0 40])
title('Baseline')


subplot(2,4,[3 4])
%plot(timev_concolidation,data.stats.activity_mean{2,1})
bar(timev_concolidation,data.stats.activity_mean{2,1},'FaceColor',[.6 .6 .6,])
hold on
plot(timev_concolidation,smooth(data.stats.activity_mean{2,1},0.4,'rloess'),'LineWidth',2);
ylim([0 40])
title('Consolidation')


subplot(2,4,5)
%plot(timev_extinction,sum(data.stats.activity_mean{3,1},3))
bar(timev_extinction,sum(data.stats.activity_mean{3,1},3),'FaceColor',[.6 .6 .6,])
hold on
plot(timev_extinction,smooth(sum(data.stats.activity_mean{3,1},3),0.4,'rloess'),'LineWidth',2);
ylim([0 40])
xlim([2 28])
title('Extinction First 10 ITI')


subplot(2,4,6)
%plot(timev_extinction,sum(data.stats.activity_mean{3,2},3),'r')
bar(timev_extinction,sum(data.stats.activity_mean{3,2},3),'FaceColor',[.6 .6 .6,])
hold on
plot(timev_extinction,smooth(sum(data.stats.activity_mean{3,2},3),0.4,'rloess'),'LineWidth',2);
ylim([0 40])
xlim([2 28])
title('Extinction Last 10 ITI')


subplot(2,4,7)
%plot(timev_extinction,sum(data.stats.activity_mean{4,1},3))
bar(timev_retrieval,sum(data.stats.activity_mean{4,1},3),'FaceColor',[.6 .6 .6,])
hold on
plot(timev_retrieval,smooth(sum(data.stats.activity_mean{4,1},3),0.4,'rloess'),'LineWidth',2);
ylim([0 40])
xlim([2 28])
title('Retrieval First 10 ITI')

subplot(2,4,8)
%plot(timev_extinction,sum(data.stats.activity_mean{4,2},3),'r')
bar(timev_retrieval,sum(data.stats.activity_mean{4,2},3),'FaceColor',[.6 .6 .6,])
hold on
plot(timev_retrieval,smooth(sum(data.stats.activity_mean{4,2},3),0.4,'rloess'),'LineWidth',2);
title('Retrieval Last 10 ITI')
ylim([0 40])
xlim([2 28])


%%

mean(data.stats.activity_mean{1,1},2)
mean(data.stats.activity_mean{2,1},2)
mean(sum(data.stats.activity_mean{3,1},3),2)
mean(sum(data.stats.activity_mean{3,2},3),2)
mean(sum(data.stats.activity_mean{4,1},3),2)
mean(sum(data.stats.activity_mean{4,2},3),2)


