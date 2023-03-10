
function sd = calculateSegNeighbors(sd, goodIntersections)

% find the segments which share the intersections of each segment

for ii = 1:length(sd)
    
    nnbs = sd(ii).intersections;
    sd(ii).segNeighbors = [];
    
    for jj = 1:length(sd)
        
        nnbs_ = sd(jj).intersections;
        sharedInts = intersect(nnbs, nnbs_);
        
        % eliminate intersections that are not on the good list
        sharedInts = intersect(sharedInts, goodIntersections);
        
        if ~isempty(sharedInts)
            
            sd(ii).segNeighbors = [sd(ii).segNeighbors, jj];
            
        end
        
    end
end