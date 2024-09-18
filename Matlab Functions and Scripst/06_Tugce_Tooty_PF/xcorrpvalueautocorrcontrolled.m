function pvalue=xcorrpvalue(series1,series2,maxi,numtests)
%xcorrpvalue Yields p value of the max value of a cross correlation
% SERIES1 and SERIES2 are the two series cross-correlated, and MAXI is the
% max of their crosscorrelogram
% NUMTESTS (optional, default 100) is the # of tests to use
% Copyright Alex Backer May 2020
numtests=100;
l1=length(series1);
count=0;
for i=1:numtests,
[xc,lags]=xcorr(series1(randperm(l1)),series2);
maxc=max(xc);
if maxc>maxi,
count=count+1;
end
end
pvalue=count/numtests;

end
