function groupArray = calcGroupArrayFromConfigArray(configArray, segsOff)

% this function lists all of the sets of connected segments from the
% configArray, each row of which lists two connected segments.
% any segments in the segsOff array are ignored.

groupArray = [];

subGroupCounter = 0;
sz = size(configArray);

go_nextSub = 1;

ind = find(configArray(:,1)>0,1,'first');

if isempty(ind)
    return;
end

while go_nextSub
    
    % find the next nonzero row in configArray
    ind = find(configArray(:,1)>0,1,'first');
    segs = configArray(ind,:);
    
    % remove segs that are off
    segs = setxor(segs, intersect(segs, segsOff));
    configArray(ind,:) = 0;
    
    if isempty(segs)
        if sum(configArray(:))==0
            go_nextSub = 0;
        end
        continue;
    end
    
    subGroupCounter = subGroupCounter+1;
    go_sub = 1;
    groupSize = 0;
    ns = length(segs);
    
    % now concatenate all of the segments that are connected to these
    % initial segments (which will form one "subgroup" and a row in
    % groupArray)
    while go_sub
        
        % add the segments to the next row in the groupArray
        groupArray(subGroupCounter,groupSize+1:groupSize+ns) = segs;
        groupSize = groupSize + ns;
        
        ind = [];
        % get the rows in which all of the current segments appear
        for ii = 1:ns
            ind = [ind, find(configArray(:)==segs(ii))];
        end
        
        if isempty(ind)
            break;
        end
        
        [row,col] = ind2sub(sz,ind);
        
        ns = length(row);
        
        % switch the column indices, since the segments connected to the
        % current segments are in the other column, whichever that is
        col_ = col;
        col(col_==1) = 2;
        col(col_==2) = 1;
        
        % find all the new segments (which are connected to the old
        % segments)
        segs = [];
        for ii = 1:ns
            segs = [segs, nonzeros(configArray(row(ii), col(ii)))];
            configArray(row(ii),:) = 0;
        end
        
        % delete off segments
        segs = setxor(segs, intersect(segs, segsOff));
        ns = length(segs);
        
        if isempty(segs)
            break;
        end
        
    end
    
    if sum(configArray(:))==0
        go_nextSub = 0;
    end
    
    
end
