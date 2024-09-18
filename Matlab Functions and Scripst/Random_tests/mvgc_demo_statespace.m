%% MVGC demo: state-space method.
%
%




%% Generate VAR test data (<mvgc_schema.html#3 |A3|>)
%
% _*Note:*_ This is where you would read in your own time series data; it should
% be assigned to the variable |X| (see below and <mvgchelp.html#4 Common
% variable names and data structures>).

%% Selected data

aa = data.lfp{7,1};

CSIT = [1 2 4];
Y = aa(:,:,CSIT);
decimate_factor = 4;

X = [];
for ii = 1:size(Y,1)
    for jj = 1:size(Y,3)
        X(ii,:,jj) = decimate(Y(ii,:,jj),decimate_factor);
    end
end

fs = parameters.decimated_srate./decimate_factor;

%% Parameters

regmode   = 'LWR';                                                            % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';                                                            % information criteria regression mode ('OLS', 'LWR' or empty for default)

morder    = 'BIC';                                                            % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = 35;                                                               % maximum model order for model order estimation

acmaxlags = [];                                                               % maximum autocovariance lags (empty for automatic calculation)

tstat     = 'F';                                                              % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;                                                             % significance level for significance test
mhtc      = 'FDRD';                                                           % multiple hypothesis test correction (see routine 'significance')

fs        = parameters.decimated_srate./mvgc.parameters.decimate_factor;      % sample rate (Hz)
fres      = 4096;                                                             % frequency resolution (empty for automatic calculation)

freqs     = sfreqs(mvgc.parameters.fres,mvgc.parameters.fs);  % Frequency vector based on frequency resolution
%% Model order estimation (<mvgc_schema.html#3 |A2|>)

% Calculate information criteria up to specified maximum model order.

ptic('\n*** tsdata_to_infocrit\n');
[AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode);
ptoc('*** tsdata_to_infocrit took ');

% Plot information criteria.

figure(1); clf;
plot_tsdata([AIC BIC]',{'AIC','BIC'},1/fs);
title('Model order estimation');

amo = 5; % actual model order

fprintf('\nbest model order (AIC) = %d\n',moAIC);
fprintf('best model order (BIC) = %d\n',moBIC);
fprintf('actual model order     = %d\n',amo);

% Select model order.

if     strcmpi(morder,'actual')
    morder = amo;
    fprintf('\nusing actual model order = %d\n',morder);
elseif strcmpi(morder,'AIC')
    morder = moAIC;
    fprintf('\nusing AIC best model order = %d\n',morder);
elseif strcmpi(morder,'BIC')
    morder = moBIC;
    fprintf('\nusing BIC best model order = %d\n',morder);
else
    fprintf('\nusing specified model order = %d\n',morder);
end

%% VAR model estimation (<mvgc_schema.html#3 |A2|>)

% Estimate VAR model of selected order from data.

ptic('\n*** tsdata_to_var... ');
[A,SIG] = tsdata_to_var(X,morder,[]);
assert(~isbad(A),'VAR estimation failed - bailing out');
ptoc;

% Report information on the estimated VAR, and check for errors.
%
% _IMPORTANT:_ We check the VAR model for stability and symmetric
% positive-definite residuals covariance matrix. _THIS CHECK SHOULD ALWAYS BE
% PERFORMED!_ - subsequent routines may fail if there are errors here. If there
% are problems with the data (e.g. non-stationarity, colinearity, etc.) there's
% also a good chance they'll show up at this point - and the diagnostics may
% supply useful information as to what went wrong.

SIG = round(SIG);

info = var_info(A,SIG);
assert(~info.error,'VAR error(s) found - bailing out');

%%
[G,info] = var_to_autocov(A,SIG);
[F_spect,fres] = autocov_to_spwcgc(G,fres,[]);

%% Granger causality calculation: time domain  (<mvgc_schema.html#3 |A13|>)

% Calculate time-domain pairwise-conditional causalities from VAR model parameters
% by state-space method [4]. The VAR model is transformed into an equivalent state-
% space model for computation. Also return p-values for specified test (F-test or
% likelihood-ratio test; this is optional - if p-values are not required, then it
% is not necessary to supply time series |X|, regression mode |regmode|, or test
% specification |tstat|).

ptic('*** var_to_pwcgc... ');
[F,pval] = var_to_pwcgc(A,SIG,X,regmode,tstat);
ptoc;

% Check for failed GC calculation

assert(~isbad(F,false),'GC calculation failed - bailing out');

% Significance-test p-values, correcting for multiple hypotheses.

sig = significance(pval,alpha,mhtc);

% Plot time-domain causal graph, p-values and significance.

figure(2); clf;
sgtitlex('Pairwise-conditional Granger causality - time domain');
subplot(1,3,1);
plot_pw(F);
title('Pairwise-conditional GC');
subplot(1,3,2);
plot_pw(pval);
title(['p-values (' tstat '-test)']);
subplot(1,3,3);
plot_pw(sig);
title(['Significant at \alpha = ' num2str(alpha)]);

%% Granger causality calculation: frequency domain  (<mvgc_schema.html#3 |A14|>)

% If not specified, we set the frequency resolution to something sensible. Warn if
% resolution is very large, as this may lead to excessively long computation times,
% and/or out-of-memory issues.

if isempty(fres)
    fres = 2^nextpow2(info.acdec); % based on autocorrelation decay; alternatively, you could try fres = 2^nextpow2(nobs);
	fprintf('\nfrequency resolution auto-calculated as %d (increments ~ %.2gHz)\n',fres,fs/2/fres);
end
if fres > 20000 % adjust to taste
	fprintf(2,'\nWARNING: large frequency resolution = %d - may cause computation time/memory usage problems\nAre you sure you wish to continue [y/n]? ',fres);
	istr = input(' ','s'); if isempty(istr) || ~strcmpi(istr,'y'); fprintf(2,'Aborting...\n'); return; end
end

% Calculate spectral pairwise-conditional causalities at given frequency
% resolution by state-space method.

ptic('\n*** var_to_spwcgc... ');
f = var_to_spwcgc(A,SIG,fres);
assert(~isbad(f,false),'spectral GC calculation failed - bailing out');
ptoc;

% Plot spectral causal graph.

figure;

sgtitlex('Pairwise-conditional Granger causality - frequency domain');
%plot_spw(F_spect,fs,[2 12]);
plot_spw(f,fs,[2 12]);

%% Granger causality calculation: frequency domain -> time-domain  (<mvgc_schema.html#3 |A15|>)

% Check that spectral causalities average (integrate) to time-domain
% causalities, as they should according to theory.

fprintf('\nfrequency-domain GC integration check... ');
Fint = smvgc_to_mvgc(f); % integrate spectral MVGCs
amax = maxabs(F+Fint)/2;
if amax < 1e-5; amax = 1; end % in case all GCs very small
mre = maxabs(F-Fint)/amax;
if mre < 1e-5
    fprintf('OK (maximum relative error ~ %.0e)\n',mre);
else
    fprintf(2,'WARNING: high maximum relative error ~ %.0e\n',mre);
end

%%
% <mvgc_demo_statespace.html back to top>
