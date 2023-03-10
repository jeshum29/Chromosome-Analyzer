function [lengths, groupsOut] = calcLengthsFromConfig(sd, segs, groups, segsOff)

% sd is segdata
% segs is a list of segments that we're considering
% groups is an array of connected segments in that group
% (each row lists a set of connected segments)

gsegs = unique(nonzeros(groups(:)));
sz = size(groups);
lengths = [];
groupsOut = [];
segsSingle = [];

for ii = 1:length(segs)
    if isempty(intersect(segs(ii), groups(:))) && isempty(intersect(segs(ii), segsOff))
        lengths(end+1) = sd(segs(ii)).length;
        groupsOut(end+1) = 0;
        segsSingle = [segsSingle, segs(ii)];
    end
end

segs = setxor(segsSingle, segs);

for ii = 1:length(segs)
    
    
    [row, col] = ind2sub(sz, find(groups==segs(ii)));
    
    if isempty(row)
        continue;
    end
    
    subsegs = unique(nonzeros(groups(row, :)));
    
    len = 0;
    for jj = 1:length(subsegs)
        len = len + sd(subsegs(jj)).length;
    end
    
    lengths(end+1) = len;
    groupsOut(end+1) = row;
    
    groups(row, :) = 0;
    
end