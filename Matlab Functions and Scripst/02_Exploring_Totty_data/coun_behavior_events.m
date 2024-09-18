a = [];
b = [];

for ii = 1:size(data.events_behavior,1)
    for jj = 1:5

        if isempty(data.events_behavior{ii,jj})
            continue
        else
            a(ii,jj) = size(data.events_behavior{ii,jj},1);
            
        end


    end
end

b(:,1) = sum(a(:,1:4),2);
