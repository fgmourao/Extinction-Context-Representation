%% Cross correlation

[acf,lags] = xcorr(data.lfp{6,3}(3,:),data.lfp{6,3}(1,:),100,'coeff');
lags = lags.*(1./parameters.decimated_srate); % convert samples to time
 
figure
plot(lags,acf); grid on;
xlabel('time lag (s)'); ylabel('correlation ({\itr})'); title('xcorr');

%% Instantaneous Frequency
fun_myfilters(data.lfp{5,1}(jj,:),parameters.decimated_srate,parameters.filter.longcutoff_1,'eegfilt',params);
instfreq(data.envelope_ALL_trials(2,1:end-1,ii),params.srate,'method','hilbert')