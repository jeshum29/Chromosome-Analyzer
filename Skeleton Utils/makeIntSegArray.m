

function intSegArray = makeIntSegArray(sd, intersections)

%  the ith row of intSegArray gives all the segments touching the
% ith intersection in the array intersections

intSegArray = [];

for ii = 1:length(intersections)
    
    counter = 0;
    
    for jj = 1:length(sd)
        
        if ~sd(jj).flag
            continue;
        end
        
        if ~isempty(intersect(intersections(ii), sd(jj).intersections))
            counter = counter + 1;
            intSegArray(ii, counter) = jj;
        end
        
    end
end