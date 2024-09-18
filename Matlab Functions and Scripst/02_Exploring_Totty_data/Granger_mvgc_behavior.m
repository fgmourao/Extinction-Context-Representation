%% Granger Causality

% 1) MVGC Multivariate Granger Causality MATLAB toolbox
% hosted at http://www.sussex.ac.uk/sackler/mvgc
% Current version is mvgc_v1.3, last updated March 2022

% L. Barnett and A. K. Seth, Granger causality for state-space models, Phys. Rev. E 91(4) Rapid Communication, 2015.
% L. Barnett and A. K. Seth, "The MVGC Multivariate Granger Causality Toolbox: A new approach to Granger-causal inference", J. Neurosci. Methods 223, pp 50-68, 2014.



% Part of this code was adapted from "mvgc_demo_statespace.m"
% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  01/2024
% Last update: 02/2024

%% Start ToolBox
startup_mvgc1 % path /Users/flavio/Documents/MATLAB/MVGC1-master
clc

clear('have_genvar_mex','mvgc_root','mvgc_version')

%% Selected data

mvgc = [];

%number of channels
mvgc.parameters.number_channels = size(data.lfp{5,1},1);

% Freezing nd Non-Freezing epochs
data_2_use = data.lfp_behavior;

% The data will be downsampled to 250 and z-scored.
mvgc.parameters.decimate_factor = 5; %
mvgc.data = [];

for ii = 1:size(data_2_use,1)
    for jj = 1:size(data_2_use,2)

        if isempty(data_2_use{ii,jj})
            continue
        end

        for cc = 1:size(data_2_use{ii,jj},1)
            for tt = 1:size(data_2_use{ii,jj},3)

                mvgc.data{ii,jj}(cc,:,tt) = zscore(decimate(data_2_use{ii,jj}(cc,:,tt),mvgc.parameters.decimate_factor));

            end
        end
    end
end


clear('data_2_use','tt','ii','jj','cc')

%% Parameters

mvgc.parameters.regmode   = 'LWR';                                                            % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
mvgc.parameters.icregmode = 'LWR';                                                            % information criteria regression mode ('OLS', 'LWR' or empty for default)

mvgc.parameters.morder    = 'BIC';                                                            % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
mvgc.parameters.momax     = 100;                                                               % maximum model order for model order estimation

mvgc.parameters.acmaxlags = [];                                                               % maximum autocovariance lags (empty for automatic calculation)

mvgc.parameters.tstat     = 'F';                                                              % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
mvgc.parameters.alpha     = 0.05;                                                             % significance level for significance test
mvgc.parameters.mhtc      = 'FDRD';                                                           % multiple hypothesis test correction (see routine 'significance')

mvgc.parameters.fs        = parameters.decimated_srate./mvgc.parameters.decimate_factor;      % sample rate (Hz)
mvgc.parameters.fres      = 4096;                                                             % frequency resolution (empty for automatic calculation)

mvgc.parameters.freqs     = sfreqs(mvgc.parameters.fres,mvgc.parameters.fs);                  % frequency vector based on frequency resolution

%% Model order estimation (<mvgc_schema.html#3 |A2|>)
% The model's order will be defined only once, based primarily on the baseline without freezing, but in the absence of this event in the first consecutive window

data_2_estimate = [];

for ii = 1:size(mvgc.data,1)
    for jj = 1:size(mvgc.data,2)

        if isempty(mvgc.data{ii,jj})
            continue

        else
            data_2_estimate = mvgc.data{ii,jj};
        end
    break
    end



end


% Calculate information criteria up to specified maximum model order.

ptic('\n*** tsdata_to_infocrit\n');
[mvgc.parameters.AIC,mvgc.parameters.BIC,mvgc.parameters.moAIC,mvgc.parameters.moBIC] = tsdata_to_infocrit(data_2_estimate,mvgc.parameters.momax,mvgc.parameters.icregmode);
ptoc('*** tsdata_to_infocrit took ');

% Plot information criteria.

mvgc.parameters.amo = 5; % actual model order

fprintf('\nbest model order (AIC) = %d\n',mvgc.parameters.moAIC);
fprintf('best model order (BIC) = %d\n',mvgc.parameters.moBIC);
fprintf('actual model order     = %d\n',mvgc.parameters.amo);

% Select model order.

if     strcmpi(mvgc.parameters.morder,'actual')
    mvgc.parameters.morder = mvgc.parameters.amo;
    fprintf('\nusing actual model order = %d\n',mvgc.parameters.morder);

    figure;
    plot_tsdata([mvgc.parameters.AIC mvgc.parameters.BIC]',{'AIC','BIC'},1/mvgc.parameters.fs);
    title(['Model order estimation - using the best actual model order = ' num2str(mvgc.parameters.morder)]);

elseif strcmpi(mvgc.parameters.morder,'AIC')
    mvgc.parameters.morder = mvgc.parameters.moAIC;
    fprintf('\nusing AIC best model order = %d\n',mvgc.parameters.morder);

    figure;
    plot_tsdata([mvgc.parameters.AIC mvgc.parameters.BIC]',{'AIC','BIC'},1/mvgc.parameters.fs);
    title(['Model order estimation - using AIC best model order = ' num2str(mvgc.parameters.morder)]);

elseif strcmpi(mvgc.parameters.morder,'BIC')
    mvgc.parameters.morder = mvgc.parameters.moBIC;
    fprintf('\nusing BIC best model order = %d\n',mvgc.parameters.morder);

    figure;
    plot_tsdata([mvgc.parameters.AIC mvgc.parameters.BIC]',{'AIC','BIC'},1/mvgc.parameters.fs);
    title(['Model order estimation - using BIC best model order = ' num2str(mvgc.parameters.morder)]);

else
    fprintf('\nusing specified model order = %d\n',mvgc.parameters.morder);

    figure;
    plot_tsdata([mvgc.parameters.AIC mvgc.parameters.BIC]',{'AIC','BIC'},1/mvgc.parameters.fs);
    title(['Model order estimation - using specified best model order = ' num2str(mvgc.parameters.morder)]);

end

clear('data_2_estimate','ii','jj');

%% Save Granger model estimation figure

% Settings
% ms = 1;
newStr = id(1:end-8);
Path        = files.FilesLoaded{1,1}(ms).folder;
name_ = strcat(Path,'/',newStr,'_Granger_model_estimation');

% save figure
saveas(gcf,name_,'png')

close all

clear('name_','newStr','path')

%% VAR model estimation (<mvgc_schema.html#3 |A2|>)

% Estimate VAR model of selected order from data.

for ii = 1:size(mvgc.data,1)
    for jj = 1:size(mvgc.data,2)

        if isempty(mvgc.data{ii,jj})
            continue
        end


        ptic('\n*** tsdata_to_var... ');
        [mvgc.parameters.A{ii,jj},mvgc.parameters.SIG{ii,jj}] = tsdata_to_var(mvgc.data{ii,jj},mvgc.parameters.morder,[]);
        assert(~isbad(mvgc.parameters.A{ii,jj}),'VAR estimation failed - bailing out');
        ptoc;

        % Report information on the estimated VAR, and check for errors.
        %
        % _IMPORTANT:_ We check the VAR model for stability and symmetric
        % positive-definite residuals covariance matrix. _THIS CHECK SHOULD ALWAYS BE
        % PERFORMED!_ - subsequent routines may fail if there are errors here. If there
        % are problems with the data (e.g. non-stationarity, colinearity, etc.) there's
        % also a good chance they'll show up at this point - and the diagnostics may
        % supply useful information as to what went wrong.

        mvgc.info{ii,jj} = var_info(mvgc.parameters.A{ii,jj},mvgc.parameters.SIG{ii,jj});
        assert(~mvgc.info{ii,jj}.error,'VAR error(s) found - bailing out');

    end
end

clear('jj','ii')

%% Granger causality calculation: frequency domain  (<mvgc_schema.html#3 |A14|>)

% If not specified, we set the frequency resolution to something sensible. Warn if
% resolution is very large, as this may lead to excessively long computation times,
% and/or out-of-memory issues.

if isempty(mvgc.parameters.fres)
    mvgc.parameters.fres = 2^nextpow2(mvgc.info.acdec); % based on autocorrelation decay; alternatively, you could try fres = 2^nextpow2(nobs);
    fprintf('\nfrequency resolution auto-calculated as %d (increments ~ %.2gHz)\n',mvgc.parameters.fres,mvgc.parameters.fs/2/mvgc.parameters.fres);
end

if mvgc.parameters.fres > 20000 % adjust to taste
    fprintf(2,'\nWARNING: large frequency resolution = %d - may cause computation time/memory usage problems\nAre you sure you wish to continue [y/n]? ',mvgc.parameters.fres);
    istr = input(' ','s'); if isempty(istr) || ~strcmpi(istr,'y'); fprintf(2,'Aborting...\n'); return; end
end

mvgc.F_spect = cell(6,5);

% Calculate spectral pairwise-conditional causalities at given frequency
% resolution by state-space method.

for ii = 1:size(mvgc.data,1)
    for jj = 1:size(mvgc.data,2)

        if isempty(mvgc.data{ii,jj})
            continue
        end

        ptic('\n*** var_to_spwcgc... ');
        mvgc.F_spect{ii,jj} = var_to_spwcgc(mvgc.parameters.A{ii,jj},mvgc.parameters.SIG{ii,jj},mvgc.parameters.fres);
        assert(~isbad(mvgc.F_spect{ii,jj},false),'spectral GC calculation failed - bailing out');
        ptoc;

    end
end

clear('ii','jj')

%% Granger causality calculation: time domain  (<mvgc_schema.html#3 |A13|>)

% Calculate time-domain pairwise-conditional causalities from VAR model parameters
% by state-space method [4]. The VAR model is transformed into an equivalent state-
% space model for computation. Also return p-values for specified test (F-test or
% likelihood-ratio test; this is optional - if p-values are not required, then it
% is not necessary to supply time series |X|, regression mode |regmode|, or test
% specification |tstat|).

for ii = 1:size(mvgc.data,1)
    for jj = 1:size(mvgc.data,2)

        if isempty(mvgc.data{ii,jj})
            continue
        end

        ptic('*** var_to_pwcgc... ');
        [mvgc.stats.F{ii,jj},mvgc.stats.pval{ii,jj}] = var_to_pwcgc(mvgc.parameters.A{ii,jj},mvgc.parameters.SIG{ii,jj},mvgc.data{ii,jj},mvgc.parameters.regmode,mvgc.parameters.tstat);
        ptoc;

    end
end

% Check for failed GC calculation
for ii = 1:size(mvgc.data,1)
    for jj = 1:size(mvgc.data,2)

        if isempty(mvgc.data{ii,jj})
            continue
        end

        assert(~isbad(mvgc.stats.F{ii,jj},false),'GC calculation failed - bailing out');

        % Significance-test p-values, correcting for multiple hypotheses.

        mvgc.stats.sig{ii,jj} = significance(mvgc.stats.pval{ii,jj},mvgc.parameters.alpha,mvgc.parameters.mhtc);

    end
end

clear('jj','ii')

%% Plot STATS -->  This code NEED TO ADAPT and UPDATE

%  -->  This code NEED TO ADAPT and UPDATE

% Plot time-domain causal graph, p-values and significance.

% ms = 1;

% % Baseline
% figure
% set(gcf, 'Position', get(0, 'Screensize'));
%
% sgtitlex('Baseline - Pairwise-conditional Granger causality - time domain');
%
% subplot(1,3,1);
% plot_pw(mvgc.stats.F{1,1});
% title('Pairwise-conditional GC');
% hb = colorbar('location','southoutside');
%
% subplot(1,3,2);
% plot_pw(mvgc.stats.pval{1,1});
% title(['p-values (' mvgc.parameters.tstat   '-test)']);
% hb = colorbar('location','southoutside');
%
% subplot(1,3,3);
% plot_pw(mvgc.stats.sig{1,1});
% title(['Significant at \alpha = ' num2str(mvgc.parameters.alpha)]);
% hb = colorbar('location','southoutside');
%
% newStr1 = id(1:end-20);
% name_1 = strcat(Path,'/',newStr1,'_Granger_Stats_baseline');
% saveas(gcf,name_1,'png') % save figure
% close all


% CS or ITI trials
% CSIT{ms} = [1 2 4]; % Choose CS and or ITI
% Session = 2; % Type 2 to CS and 3 to ITI
% sub_idx = reshape((1:length(CSIT{ms})*3),[],3);
%
% figure
% set(gcf, 'Position', get(0, 'Screensize'));
%
% sgtitlex('Pairwise-conditional Granger causality - time domain - CS Trials');
%
%
% for ii = 1:size(sub_idx,1)
%     s = subplot(3,size(sub_idx,1),sub_idx(ii,1));
%     plot_pw(mvgc.stats.F{Session,CSIT{ms}(ii)});
%     title({['F-test'],['Trial: ' num2str(CSIT{ms}(ii))]});
%     %clim([0 1])
%
%     if ii == size(sub_idx,1)
%         s2Pos = get(s,'position');
%         hb = colorbar('location','eastoutside');
%         set(s,'position',s2Pos);
%     end
%
% end
%
% for ii = 1:size(sub_idx,1)
%     s = subplot(3,size(sub_idx,1),sub_idx(ii,2));
%     plot_pw(mvgc.stats.pval{Session,CSIT{ms}(ii)});
%     title(['p-values (' mvgc.parameters.tstat   '-test)']);
%
%     if ii == size(sub_idx,1)
%         s2Pos = get(s,'position');
%         hb = colorbar('location','eastoutside');
%         set(s,'position',s2Pos);
%     end
%
% end
%
% for ii = 1:size(sub_idx,1)
%     s = subplot(3,size(sub_idx,1),sub_idx(ii,3));
%     plot_pw(mvgc.stats.sig{Session,CSIT{ms}(ii)});
%     title(['Significant at \alpha = ' num2str(mvgc.parameters.alpha)]);
%
%     if ii == size(sub_idx,1)
%         colorbar
%         s2Pos = get(s,'position');
%         hb = colorbar('location','eastoutside');
%         set(s,'position',s2Pos);
%     end
%
% end
%
% % Save
% newStr1 = id(1:end-8);
% Path    = files.FilesLoaded{1,1}(ms).folder;
% name_1  = strcat(Path,'/',newStr1,'_Granger_Stats_CSTrials');
%
% saveas(gcf,name_1,'png') % save figure
%
% close all
%
% clear('ii','s','s2Pos','hb','sub_idx','name_1','newStr1')

%% Granger causality calculation: frequency domain -> time-domain  (<mvgc_schema.html#3 |A15|>)

% -->  This code NEED TO ADAPT and UPDATE

% % Check that spectral causalities average (integrate) to time-domain
% % causalities, as they should according to theory.
% 
% % Define band frequency
% % 3 - 6 Hertz
% steps                          = diff(mvgc.parameters.freqs); % according to the fft time window
% mvgc.parameters.frex_3_6Hz     = 3:steps(1):6;
% mvgc.parameters.frex_idx_3_6Hz = dsearchn(mvgc.parameters.freqs,mvgc.parameters.frex_3_6Hz');
% 
% for tt = 1:size(mvgc.data,1)
% 
%     if isempty(mvgc.data{tt,1})
%         continue
%     end
% 
%     for ii = 1:size(mvgc.data{tt,1},3)
% 
% 
%         fprintf('\nfrequency-domain GC integration check... ');
%         mvgc.Fint_3_6Hz{tt,ii} = smvgc_to_mvgc(mvgc.F_spect{tt,ii}(:,:,mvgc.parameters.frex_idx_3_6Hz)); % integrate spectral MVGCs
% 
% 
%         mvgc.parameters.amax = maxabs(mvgc.stats.F{tt,ii}+mvgc.Fint_3_6Hz{tt,ii})/2;
% 
%         if mvgc.parameters.amax < 1e-5; mvgc.parameters.amax = 1; end % in case all GCs very small
%         mvgc.parameters.mre = maxabs(mvgc.stats.F{tt,ii}-mvgc.Fint_3_6Hz{tt,ii})/mvgc.parameters.amax;
%         if mvgc.parameters.mre < 1e-5
%             fprintf('OK (maximum relative error ~ %.0e)\n',mvgc.parameters.mre);
%         else
%             fprintf(2,'WARNING: high maximum relative error ~ %.0e\n',mvgc.parameters.mre);
%         end
% 
%     end
% end
% 
% 
% % 6 - 9 Hertz
% steps                          = diff(mvgc.parameters.freqs); % according to the fft time window
% mvgc.parameters.frex_6_9Hz     = 6:steps(1):9;
% mvgc.parameters.frex_idx_6_9Hz = dsearchn(mvgc.parameters.freqs,mvgc.parameters.frex_6_9Hz');
% 
% for tt = 1:size(mvgc.data,1)
% 
%     if isempty(mvgc.data{tt,1})
%         continue
%     end
% 
%     for ii = 1:size(mvgc.data{tt,1},3)
% 
%         fprintf('\nfrequency-domain GC integration check... ');
%         mvgc.Fint_6_9Hz{tt,ii} = smvgc_to_mvgc(mvgc.F_spect{tt,ii}(:,:,mvgc.parameters.frex_idx_6_9Hz)); % integrate spectral MVGCs
% 
% 
%         mvgc.parameters.amax = maxabs(mvgc.stats.F{tt,ii}+mvgc.Fint_6_9Hz{tt,ii})/2;
% 
%         if mvgc.parameters.amax < 1e-5; mvgc.parameters.amax = 1; end % in case all GCs very small
%         mvgc.parameters.mre = maxabs(mvgc.stats.F{tt,ii}-mvgc.Fint_6_9Hz{tt,ii})/mvgc.parameters.amax;
%         if mvgc.parameters.mre < 1e-5
%             fprintf('OK (maximum relative error ~ %.0e)\n',mvgc.parameters.mre);
%         else
%             fprintf(2,'WARNING: high maximum relative error ~ %.0e\n',mvgc.parameters.mre);
%         end
% 
%     end
% end
% 
% clear('tt','ii','steps')

%% Save data

% Settings
%ms = 1;
newStr = id(1:end-8);
%Path    = files.FilesLoaded{1,1}(ms).folder;
Path = '/Users/flavio/Desktop';
%name_1  = strcat(Path,'/',newStr1,'_Granger_Spec_Mean_CS_Trials_2');
name  = strcat(Path,'/',newStr,'_Granger_Freezing_NonFreezing');


% Save data
save(name,'mvgc','-v7.3')

clear('name','newStr','path')


%% last update 25/02/2024
%  listening:
