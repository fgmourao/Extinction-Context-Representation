
filtt = [];
temp = [];
filtt.parameters.filters = [3.5 5;7 9];
% Baseline
for ii = 1:size(data.lfp{5,1},1)
    for jj = 1:size(filtt.parameters.filters,1)
        temp{1,jj}(ii,:) = eegfilt(data.lfp{5,1}(ii,:),parameters.decimated_srate,filtt.parameters.filters(jj,1),filtt.parameters.filters(jj,2)); % just filtering
        filtt.amplitude{1,jj}(ii,:) = abs(hilbert(temp{1,jj}(ii,:))); % getting the amplitude envelope
        filtt.phase{1,jj}(ii,:) = angle(hilbert(temp{1,jj}(ii,:))); % getting the angles
    end
end

%%
tt = 1;
tbl = 2000;
ev = 4;

timev = linspace(0,14,length(data.events{ev, 1}(tt,1)-tbl:data.events{ev, 1}(tt,2)+tbl));
figure
subplot(2,1,1)
plot(timev,temp{1, 2}(3,data.events{ev, 1}(tt,1)-tbl:data.events{ev, 1}(tt,2)+tbl))
hold
plot(timev,temp{1, 2}(1,data.events{ev, 1}(tt,1)-tbl:data.events{ev, 1}(tt,2)+tbl))
subplot(2,1,2)
plot(timev,temp{1, 2}(3,data.events{ev, 1}(tt,1)-tbl:data.events{ev, 1}(tt,2)+tbl))
hold
plot(timev,temp{1, 2}(2,data.events{ev, 1}(tt,1)-tbl:data.events{ev, 1}(tt,2)+tbl))
