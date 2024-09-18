% 
% Since the fiber photometry data has a sampling frequency (20 Hz) different 
% from the video freeze data (movement behavior; 30 Hz), the data were interpolated here. 
% The fiber photometry data were normalized to deltaf/f using PMAT and exported as *.csv.


for ii = 2:size(Exposurerawactivity,2)

    xq = linspace(PhotometryDeltaFtraces(1,1),PhotometryDeltaFtraces(end,1),length(Exposurerawactivity));
    data_iterp_(ii,:) = interp1(PhotometryDeltaFtraces(:,1),PhotometryDeltaFtraces(:,ii),xq)';

end

% Save as a table
variables2save = {'Time(s)','F/f','Movement'};
spreadsheet1 = array2table(SP6_,'VariableNames',variables2save);

% Save *.csv
newStr = 'SP6';
path = '/Users/flavio/Desktop';
name1 = strcat(path,'/',newStr,'_data_extinction','.csv');

writetable(spreadsheet1,name1);
