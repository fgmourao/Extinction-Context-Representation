%gauss_filt = fspecial('gaussian',[1 400],1);


cell_z = [];

cell_z{1,1} = FearNondecrementcells(:,1:400);
cell_z{2,1} = FearNondecrementcells(:,402:801);

cell_z{3,1} = Fearcells(:,1:400);
cell_z{4,1} = Fearcells(:,402:802);

cell_z{5,1} = Extinctioncells(:,1:400);
cell_z{6,1} = Extinctioncells(:,402:801);


%%
gauss_filt = fspecial('gaussian',[1 400],1);
cell_z_filter_norm = [];

for jj = 1:size(cell_z,1)
    for ii = 1:size(cell_z{jj,1},1)
        cell_z_filter_norm{jj,1}(ii,:) = conv(mapminmax(cell_z{jj,1}(ii,:),0,1),gauss_filt,'same');
    end
   
end

cell_z_filter = [];

for jj = 1:size(cell_z,1)
    for ii = 1:size(cell_z{jj,1},1)
        cell_z_filter{jj,1}(ii,:) = cell_z{jj,1}(ii,:);
    end
end

%%

cell_z_binary = [];

for jj = 1:size(cell_z,1)
    cell_z_binary{jj,1} = zeros(size(cell_z{jj,1}));
end

for jj = 1:size(cell_z,1)
    for ii = 1:size(cell_z{jj,1},1)
        temp = find(cell_z{jj,1}(ii,:)>1.93);
        cell_z_binary{jj,1}(ii,temp) = 1;
    end
end



%% cell_z_filter_norm_sum

cell_z_sum = [];

for jj = 1:size(cell_z,1)
    cell_z_sum{jj,1} = sum(cell_z_binary{jj,1}(:,102:106),2);

end

%     cell_z_filter_norm_sum{jj,2} = std(cell_z_filter_norm{jj,1},[],2,'omitnan')./size(cell_z_filter_norm{jj,1},1);
% 
%     curve1_ext{ii,1} = data_2_plot_mscohr_retrieval_mean_SEM{ii, 1} + data_2_plot_mscohr_retrieval_mean_SEM{ii, 2};
%     curve2_ext{ii,1} = data_2_plot_mscohr_retrieval_mean_SEM{ii, 1} - data_2_plot_mscohr_retrieval_mean_SEM{ii, 2};
%     data_2_plot_mscohr_retrieval_mean_SEM{ii,3} = cat(2,curve1_ext{ii,1},flip(curve2_ext{ii,1},2));
%     data_2_plot_mscohr_retrieval_mean_SEM{ii,4} = [freq_v_mscohr', fliplr(freq_v_mscohr')];
% 
% 
% 
% end

temp = [cell_z_sum{1,1} cell_z_binary{1,1}];
cell_z_sorted{1,1} =  sortrows(temp,1,'descend');

temp = [cell_z_sum{1,1} cell_z_binary{2,1}];
cell_z_sorted{2,1} =  sortrows(temp,1,'descend');

temp = [cell_z_sum{3,1} cell_z_binary{3,1}];
cell_z_sorted{3,1} =  sortrows(temp,1,'descend');

temp= [cell_z_sum{3,1} cell_z_binary{4,1}];
cell_z_sorted{4,1} =  sortrows(temp,1,'descend');

temp= [cell_z_sum{6,1} cell_z_binary{5,1}];
cell_z_sorted{5,1} =  sortrows(temp,1,'descend');

temp = [cell_z_sum{6,1} cell_z_binary{6,1}];
cell_z_sorted{6,1} =  sortrows(temp,1,'descend');


%%
% early 10 --> 1 - 400 samples
% late  10 --> 402 - 802 samples

%data2plot = cell_z;
data2plot = cell_z_sorted;

time_v = linspace(-5,15,400);

figure
set(gcf,'color','white')
box 'off'
hold on

t = {'FearNondecrementcells early','FearNondecrementcells late','Fearcells early','Fearcells late','Extinctioncells early','Extinctioncells late'};

for ii = 1:size(data2plot,1)

    subplot(3,2,ii)
    imagesc(time_v,1:size(data2plot{ii,1},1),data2plot{ii,1}(:,2:end));
    %contourf(time_v,1:size(data2plot{ii,1},1),data2plot{ii,1}(:,1:400),1240,'linecolor','none');
    title(t{1,ii});
    a = gca;
    a.YColor = 'w';
    a.YTick = [];
    xlim([-2 12])
    clim([0 1])
    colorbar
    xlabel('\fontsize{12}Time (s)');
    %ylabel('\fontsize{12}Cells (1 -> 57)');

    %c = colorbar;
    colormap parula
    colorbar off
    %set(get(c,'ylabel'),'string','\fontsize{12} Zscore','Rotation',270);
    %set(c,'XTickLabel',{'0',' ',' ',' ',' ','5'});
    %view(0,90)
%     py_path = "~/anaconda3/envs/Python_3_10/bin/python";
%     Py_map = getPyPlot_cMap('viridis', [], [], py_path);
%     colormap(Py_map)


end


%%

%data2plot = cell_z;
data2plot = cell_z_sorted;

time_v = linspace(-5,15,400);

figure
set(gcf,'color','white')
sc = [1,1,960,1200];
set(gcf, 'Position', sc);
hold on

for ii = 1:size(data2plot,1)

    subplot(3,2,ii)
    bar(time_v,sum(data2plot{ii,1}(:,2:401),1),'FaceColor',[.6 .6 .6,])
%     hold on
%     plot(time_v,smooth(sum(data2plot{ii,1}(:,2:401),1),0.013,'rloess'),'LineWidth',2);
    xlim([9.8 13])
    ylim([0 14])
    xline(0,'--k', 'LineWidth',1)
    xline(10,'--k','LineWidth',1)
    xlabel('\fontsize{12}Time (s)');
    ylabel('\fontsize{12}Bin Count');
    box 'off'


    title(t{1,ii});


end



