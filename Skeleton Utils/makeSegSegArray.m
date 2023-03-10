

function [segSegArray, intersections] = makeSegSegArray(sd, segs, goodIntersections)

% the ith row of segSegArray lists the segments that neighbor the ith
% segment in the segs array

% intersections is a list of all the good intersections that the segments in the
% segs array touch

intersections = [];
segSegArray = [];

for ii = 1:length(segs)
    
    if ~sd(segs(ii)).flag
        segSegArray(ii,1) = NaN;
        continue;
    end
    
    intersections = [intersections, sd(segs(ii)).intersections'];
    counter = 0;
    
    for jj = 1:length(sd(segs(ii)).segNeighbors);
        if sd(sd(segs(ii)).segNeighbors(jj)).flag
            counter = counter+1;
            segSegArray(ii, counter) = sd(segs(ii)).segNeighbors(jj);
        end
    end
    
    
end

intersections = intersect(unique(intersections), goodIntersections);


