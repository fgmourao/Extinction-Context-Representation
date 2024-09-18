Rat1Cells = Rat1Cells'; 
timev = Rat1Cells(1,:);
Rat1Cells(1,:) = [];

x = Rat1Cells;

[x_normalized, PS]=mapminmax(x,0,1);
% norm_data = (bla - minVal) / ( maxVal - minVal )
% your_original_data = minVal + norm_data.*(maxVal - minVal)

%x_normalized_z = zscore(x_normalized,0,1);
%%
figure
set(gcf,'color','white')
box 'off'
hold on
% Choose channels to plot
channels = 1:57;

% factor
factor = (channels)'*1;
r = plot(1:length(x_normalized(channels,:)), bsxfun(@plus, x_normalized(channels,:), factor),'Color','[0.3, 0.3, 0.3]','linew',1);
a = gca;
a.YColor = 'w';
a.YTick = [];
a.XLim = [0  timev(end)];
box off
%xline(180)
%xlim([0 300])


%%
test = zscore(x_normalized);

neuronts  = fun_myfilters(Exp3Day2,5,[1 0.5],'eegfilt',[]);

neuronts = zscore(bsxfun(@rdivide,Exp3Day2,mean(Exp3Day2,2)));
neuronts = zscore(Exp3Day2);



contourf(timev(1:300),1:57,(x_normalized(:,1:300)),80,'linecolor','none');
colorbar
caxis([-1 1])
