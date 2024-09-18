%%

data_f_iterp_ = cell(size(data_f));

for jj = 1:size(data_f,1)
    for ii = 1:size(data_f,2)

        xq = linspace(time_F(1,1),time_F(1,end),length(behavior_epochs_mean{jj, ii}));
        data_f_iterp_{jj,ii} = interp1(time_F,data_f{jj,ii},xq);

    end
end

%%
acf = cell(size(data_f));
lags = cell(size(data_f));
p = cell(size(data_f));

for jj = 1:size(data_f,1)
    for ii = 1:size(data_f,2)

        [acf{jj,ii},lags{jj,ii}] = xcorr(zscore(data_f_iterp_{jj,ii}),zscore(behavior_epochs_mean{jj, ii}),[],'coeff');
        lags{jj,ii} = lags{jj,ii}.*(1./30); % convert samples to time

    end
end


numshf  = 1000; % number of shuffled segments
nsurrog = 1000; % number of rearrangements

data_shuf{jj,ii} = cell(size(data_f));
acf_shuffle = [];

for ss = 1:nsurrog
    for jj = 1:size(data_f,1)
        for ii = 1:size(data_f,2)

            data_shuf{jj,ii} = shuffle_esc(behavior_epochs_mean{jj, ii},30,numshf);

        end
    end

    for jj = 1:size(data_f,1)
        for ii = 1:size(data_f,2)

            [acf_shuffle{jj,ii}(:,:,ss),~] = zscore(xcorr(zscore(data_f_iterp_{jj,ii}),zscore(data_shuf{jj,ii}),[],'coeff'));

        end
    end


end

z_observed = [];
p_value = [];

for jj = 1:size(data_f,1)
    for ii = 1:size(data_f,2)

        z_observed{jj,ii} = max(acf{jj,ii})...
            - mean(squeeze(mean(acf_shuffle{jj,ii}(1,61:81,:))))./...
            std(squeeze(mean(acf_shuffle{jj,ii}(1,61:81,:))),[],1); % z real(observed) value


        p_value{jj,ii} = 0.5 * erfc(z_observed{jj,ii}  ./ sqrt(2)); % p value. Similar to Matlab function: normcdf(-z) two-tailed

    end
end



%% Plots


figure
set(gcf,'color','w');

for cc = 1:size(acf,2)
    subplot(2,8,cc)
    hold on
    %plot(lags{1,cc},acf{1,cc},'color',[.6 .6 .6],'linew',2); grid on;
    plot(lags{1,cc},acf{1,cc},'color',[.8 0 0 .5],'linew',2); grid on;
    xlabel('time lag (s)'); ylabel('correlation ({\itr})');
    title(['SP' num2str(cc)]);
    %xlim([-2 2])
    ylim([-.5 1])

    subplot(2,8,cc+8)
    h = histogram(squeeze(mean(acf_shuffle{1,cc}(1,61:81,:),2)),50);
    h.FaceColor = [1 1 1];
    h.EdgeColor = 'k';
    xlim([-4 4])
    ylim([0 60])
    
    hold on
    stem(z_observed{1,cc},20,'linew',2);
    xlabel(['p = ' num2str(p_value{1,cc})])


end

subplot(2,8,8)
%boundedline(lags{1,1},mean(cat(3,acf{1,:}),3),(std(cat(3,acf{1,:}),[],3))./size(cat(3,acf{1,:}),3),'k'), grid on; % baseline
hold on
boundedline(lags{1,1},mean(cat(3,acf{1,:}),3),(std(cat(3,acf{1,:}),[],3))./size(cat(3,acf{1,:}),3),'cmap', hot(4)), grid on; % baseline

xlabel('time lag (s)'); ylabel('correlation ({\itr})'); title('xcorr mean');
%xlim([-2 2])
ylim([-.5 1])


%%

series1 = zscore(data_f_iterp_{1,1});
series2 = zscore(behavior_epochs_mean{1, 1});
p=xcorrpvalue(series1,series2,max(abs(xcorr(series1,series2,'coeff'))))



