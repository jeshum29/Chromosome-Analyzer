function flag = searchForLoops(configGroupArray, intSegArray)

for ii = 1:size(configGroupArray,1)
    
    flag(ii) = 0;
    segs = configGroupArray(ii,:);
    
    if length(segs)==1
        continue;
    end
            
    % look for an intersection that's shared
    % by three segments. This corresponds to 
    % a loop intersected by a segment
    used_ints = [];
    
    for jj = 1:size(intSegArray,1)
        if length(intersect(intSegArray(jj,:), segs)) == 3
            flag(ii) = 1;
        end
    end
    
    if flag(ii)
        continue;
    end
    
    
    %check for a simple loop formed by two segs
    if length(segs)==2
        inds_1 = find(intSegArray==segs(1));
        inds_2 = find(intSegArray==segs(2));
        if length(inds_1)==2 && length(inds_2)==2
            if sum(inds_1==inds_2)==2
                flag(ii) = 1;
                continue;
            end
        end
    end

    
    
    % look for a simple loop formed by three or more segments
    
    seg = segs(1);
    int = [];
    for jj = 1:length(segs)
        
        [int_new,cc] = ind2sub(size(intSegArray), find(intSegArray==seg));
        
        if isempty(int)
            int = int_new(1);
        else
            int = setxor(int_new, int);
        end
        
        next_seg = setxor(intSegArray(int,:), seg);
        seg = intersect(next_seg, segs);
        
    end
    
    % if we ended up where we started
    if seg==segs(1)
        flag(ii) = 1;
        continue;
    end

    
end

    
    
    
