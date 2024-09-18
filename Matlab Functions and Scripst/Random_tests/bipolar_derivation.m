%% Bipolar derivation (optinal)


%                           - CHANNELS MAP - Original monopolar -

% Extinction
% Animal 7_1
% mPFC PL -> 01 - 6
% mPFC PL -> 07 - 12
% HPC     -> 13 - 22

%                           - CHANNELS MAP - Bipolar -

% mPFC PL -> 1:5;
% mPFC PL -> 6:10;
% HPC     -> 11:19;

%%

for ii = 1:size(data.lfp,1)
    for jj = 1:size(data.lfp,2)

        if isempty(data.lfp{ii,jj}) % The raw data has never been filtered until now.
            continue
        end

        data.lfp{ii,jj} = diff(data.lfp{ii,jj},1,1);
        data.lfp{ii,jj}([6 12],:,:) = [];
    end
end














