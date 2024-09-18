timef_extinction = [];
for ii = 1:size(corr_extic.time_freezing_CS,1)
    temp = corr_extic.time_freezing_CS{ii,1};
    timef_extinction = [timef_extinction temp];
end

timef_retrieval = [];
for ii = 1:size(corr_retrieval.time_freezing_CS,1)
    temp = corr_retrieval.time_freezing_CS{ii,1};
    timef_retrieval = [timef_retrieval temp];
end

figure
set(gcf,'color','w');
sc = [1,1,960,1200];
set(gcf, 'Position', sc);

hold on
b1 = histogram(timef_extinction,500,'FaceAlpha',.6);

b1(1).FaceColor = [.6, 0, 0];
b1(1).EdgeColor = [.6, 0, 0];
b1.BinLimits = [0 50];

figure
set(gcf,'color','w');
sc = [1,1,960,1200];
set(gcf, 'Position', sc);

b2 = histogram(timef_retrieval,35,'FaceAlpha',.6);
b2.BinEdges = b1.BinEdges;
b2(1).FaceColor = [.2, .2, 1];
b2(1).EdgeColor = [.2, .2, 1];
b2.BinLimits = [0 50];
xlabel('\fontsize{11}Time (s)');
ylabel('\fontsize{11}Frequency');
