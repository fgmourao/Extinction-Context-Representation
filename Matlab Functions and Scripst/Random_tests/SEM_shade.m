a = squeeze(data_2_plot_extinction_mean_SEM{3, 1}(2,3,:))';
s = squeeze(data_2_plot_extinction_mean_SEM{3, 2}(2,3,:))';

curve1_ext = a + s;
curve2_ext = a - s;
x2 = [freq_v', fliplr(freq_v')];
b = [curve1_ext, fliplr(curve2_ext)];

fill(x2, b, 'k','FaceAlpha',0.6,'EdgeColor','none');
hold on
plot(freq_v,a,'linew', 2,'Color','k')