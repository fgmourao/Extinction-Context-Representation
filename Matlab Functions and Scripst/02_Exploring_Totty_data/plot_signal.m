hilb_extinction = [];
temp_extinction = [];

hilb_extinction.parameters.filters = [3 6;6 8;7 10];

% Baseline
for ii = 1:size(data.lfp{6,1},1)
    for jj = 1:size(hilb_extinction.parameters.filters,1)

        temp_extinction{1,jj}(ii,:) = eegfilt(data.lfp{6,1}(ii,10000:170000),parameters.decimated_srate,hilb_extinction.parameters.filters(jj,1),hilb_extinction.parameters.filters(jj,2)); % just filtering
        hilb_extinction.amplitude{1,jj}(ii,:) = abs(hilbert(temp_extinction{1,jj}(ii,:))); % getting the amplitude envelope
        hilb_extinction.phase{1,jj}(ii,:) = angle(hilbert(temp_extinction{1,jj}(ii,:))); % getting the angles

    end

end


%%

hilb_retrieval = [];
temp_retrieval = [];

hilb_retrieval.parameters.filters = [3 6;6 8;7 10];

% Baseline
for ii = 1:size(data.lfp{6,1},1)
    for jj = 1:size(hilb_retrieval.parameters.filters,1)

        temp_retrieval{1,jj}(ii,:) = eegfilt(data.lfp{6,1}(ii,10000:170000),parameters.decimated_srate,hilb_retrieval.parameters.filters(jj,1),hilb_retrieval.parameters.filters(jj,2)); % just filtering
        hilb_retrieval.amplitude{1,jj}(ii,:) = abs(hilbert(temp_retrieval{1,jj}(ii,:))); % getting the amplitude envelope
        hilb_retrieval.phase{1,jj}(ii,:) = angle(hilbert(temp_retrieval{1,jj}(ii,:))); % getting the angles

    end

end