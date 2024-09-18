% Check channels

figure
set(gcf,'color','w');
set(gcf, 'Position', get(0, 'Screensize'));

%sgtitle('Exposure CxB')
for ii = 1:size(data.lfp{1,1},1)
    subplot(8,4,ii)
    plot(data.lfp{1,1}(ii,:))
    ylim([-600 600])
    title(['channel: ' num2str(ii)])

end


%% Save

% Settings
% ms = 1;
newStr1 = id(1:end-8);
Path        = files.FilesLoaded{1,1}(ms).folder;
name_1 = strcat(Path,'/',newStr1,'_check_channels');

% save figure
saveas(gcf,name_1,'png')

close all

clear('name_1','newStr1','path') 
