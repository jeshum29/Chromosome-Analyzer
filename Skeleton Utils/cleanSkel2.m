function [sk_clean, flag] = cleanSkel2(sk)

%----------------------------------------------------------
%
% function to optimize a "noisy" skeleton by maximizing length
% and minimizing curvature
%
% Keith Cheveralls
% Dernburg Lab
% December 2011
%
%----------------------------------------------------------

sk_clean = [];
flag = [];

if sum(sk(:)) < 3
    sk_clean = sk;
    'cleanskel2 warning: very short skeleton';
    return;
end

skL_tmp = bwlabeln(sk, 26);

if max(skL_tmp(:)) > 1
    'cleanskel2 warning: skeleton is discontiguous';
    return;
end

sk = deleteSpur(sk,1);

[sd, skL, intL] = make_skel_data(sk);

if length(sd)==1
    sk_clean = sk;
    return;
else
    sk_clean = 0*sk;
end

% try to reduce number of segments if there are way too many

numIntersections = max(intL(:));
flag = 0;
if numIntersections > 6
    sk = deleteSpur(sk,2);
    [sd, skL, intL] = make_skel_data(sk);
    numIntersections = max(intL(:));
    
    if numIntersections > 6
        sk = deleteSpur(sk,3);
        [sd, skL, intL] = make_skel_data(sk);
        numIntersections = max(intL(:));
        
        if numIntersections > 6
            sk = deleteSpur(sk,4);
            [sd, skL, intL] = make_skel_data(sk);
            numIntersections = max(intL(:));
            
            if numIntersections > 6
                'cleanSkel2 warning: too many intersections!';
                flag = 1;
            end
        end
    end
end

%----------------------------------------------------------
%
% begin by determining segment neighbors from intersection
% and also seg-seg and int-seg interaction arrays
%
%----------------------------------------------------------

goodIntersections = 1:numIntersections;
allSegs = 1:length(sd);

%----------------------------------------------------------
%
% these three lines calculate all of the seg-seg and int-seg
% interactions using the list of good intersections (i.e., ones
% we haven't turned off) and using the sd.flag field to ignore
% "off" segments.
%
% first, update the segNeighbors field in sd to reflect deleted
% intersections
%
%----------------------------------------------------------

sd = calculateSegNeighbors(sd, goodIntersections);

%----------------------------------------------------------
%
% from the segNeighbors field, make a seg-seg interaction array.
% the ith row lists the segments that neighbor the ith segment.
% if the segment is off, the ith row is NaNs.
%
% intersections is a list of all the intersections that are used to
% make segSegArray. The intersections not appearing in
% goodIntersections are implicitly absent since they are excluded
% earlier in calculateSegNeighbors.
%
%----------------------------------------------------------

[segSegArray, intersections] = makeSegSegArray(sd, allSegs, goodIntersections);

%----------------------------------------------------------
%
% make an int-seg array. the ith row gives the segments that
% neighbor the intersection given by intersections(ii).
%
%----------------------------------------------------------

intSegArray = makeIntSegArray(sd, intersections);


%----------------------------------------------------------
%
% check for intersections between more than three segments.
% for now, ignore these by turning off the shortest otherwise
% unconnected of the four segments
%
%----------------------------------------------------------

if size(intSegArray,2)>=4
    
    count = sum(~~intSegArray, 2);
    inds = find(count>=4);
    
    for ii = 1:length(inds)
        
        segs = intSegArray(inds(ii),:);
        isolatedSeg = [];
        len = [];
        
        for jj = 1:length(segs)
            
            % if the four-way is the segment's only intersection
            if length(sd(segs(jj)).intersections)==1
                
                isolatedSeg = [isolatedSeg, segs(jj)];
                len = [len, sd(segs(jj)).length];
                
            end
            
        end
        
        if isempty(isolatedSeg)
            % turn off the intersection. this is a drastic step
            goodIntersections = setdiff(goodIntersections, inds(ii));
            'warning: a nontrivial four-way intersection was detected';
            continue;
        else
            [mm, jj_ind] = min(len);
        end
        
        if ~sd(isolatedSeg(jj_ind)).hardFlag
            % turn off the shortest dangling segment
            sd(isolatedSeg(jj_ind)).flag = 0;
        else
            % turn off the intersection. this is a drastic step
            goodIntersections = setdiff(goodIntersections, inds(ii));
            'warning: a nontrivial four-way intersection was detected';
        end
        
    end
    
    % recalculate the global interaction arrays
    sd = calculateSegNeighbors(sd, goodIntersections);
    [segSegArray, intersections] = makeSegSegArray(sd, allSegs, goodIntersections);
    intSegArray = makeIntSegArray(sd, intersections);
    
end



%----------------------------------------------------------
%
% there might be some intersections that have only two segments,
% because we turn off bad (short) segments
%

real_intersections = sum(~~intSegArray, 2)==3;

threeWayIntersections = intSegArray(real_intersections, :);
twoWayIntersections = intSegArray(~real_intersections, :);

twoWayTmp = [];
for jj = 1:size(twoWayIntersections,1)
    twoWayTmp(jj, 1:2) = nonzeros(twoWayIntersections(jj,:));
end
twoWayIntersections = twoWayTmp;

%
%----------------------------------------------------------

%----------------------------------------------------------
%
% now explore all possible sets of two-way intersections that
% can be constructed from the three-ways
%
%----------------------------------------------------------

N = sum(real_intersections);
numPerms = 4^N;
clear numSegsOff;
clear masterLengths;

if N > 5
    
    'TOO MANY INTERSECTIONS'
    
    sk_clean = [];
    flag = [];
    return;
    
    
end

for nn = 1:(numPerms+1)
    
    if ~mod(nn,1000)
        nn;
    end
    
    hardFail = 0;
    
    % use the digits of a base-four number to choose which two of three
    % segments at each intersection are "on" (fourth option: none)

    try
    [configArray, unconnectedSegs] = makeConfigArray(threeWayIntersections, nn-1, N);
    catch
        '';
    end
    
    if hardFail
        continue;
    end
    
    if ~isempty(twoWayIntersections)
        configArray = cat(1, twoWayIntersections, configArray);
    end
    
    segsOff = [];
        
    %config group array lists all the groups of connected segments
    configGroupArray = calcGroupArrayFromConfigArray(configArray, segsOff);
    
%     loopFlag = searchForLoops(configGroupArray, intSegArray);

try
    [segLengths, segGroups] = calcLengthsFromConfig(sd, allSegs, configGroupArray, segsOff);
    catch
        '';
    end
    
    allSegLengths(nn, 1:length(segLengths)) = segLengths;
    allSegGroups(nn,1:length(segGroups)) = segGroups;

    
end

% choose the configuration that has the longest segment that isn't a loop
found_non_loop = 0;
[mmax, ind] = max(allSegLengths(:));

for nn = 1:size(allSegLengths,1)
    
    flag = 0;
    [rr,cc] = ind2sub(size(allSegLengths), ind);

    [configArray, unconnectedSegs] = makeConfigArray(threeWayIntersections, rr-1, N);
    configGroupArray = calcGroupArrayFromConfigArray(configArray, segsOff);
    good_segs = configGroupArray(allSegGroups(rr,cc), :);

    flag = loopFlagger(good_segs, intSegArray);

    % if there's a loop, check the next longest configuration
    if flag
        allSegLengths(rr,cc) = 0;
        [mmax, ind] = max(allSegLengths(:));
    else
        found_non_loop = 1;
        break;
    end

end

if ~found_non_loop
    'warning: no non-loop configs found!';
end

if 0
    vizSkel_Config(skL, intL, sd, allSegs, threeWayIntersections, twoWayIntersections, rr-1, N)
end

% recover the segments that were connected together to make this longest
% segment


for ii = 1:length(good_segs)
    
    % turn on each good seg
    sk_clean(skL==good_segs(ii)) = 1;
    
    % find the intersections each seg touches
    ind = find(intSegArray(:)==good_segs(ii));
    [rr,cc] = ind2sub(size(intSegArray), ind);
    

    % if the segment is part of a loop, there will be more than one
    % intersection
    for jj = 1:length(rr) 
        % turn on the intersections
        sk_clean(intL==rr(jj)) = 1;
    end
    
end


sk_clean = skel3(double(sk_clean));

sk_clean = preprocess_skel(sk_clean);

'';





function sk = preprocess_skel(sk)

% remove small spurs that we know are bad

[sd, skL, intL] = make_skel_data(sk);

for ii = 1:length(sd)
    
    if sd(ii).length == 1
        sk(skL==ii) = 0;
    end
    
end

sk = skel3(double(sk));





function [skL, intL] = label_skel(sk)

int = find_intersections_3(double(sk), 26);
skL = sk;
skL(int>2) = 0;
skL = bwlabeln(skL,26);
intL = bwlabeln(int>2,26);

function [sd, skL, intL] = make_skel_data(sk)

sz = size(sk);
[skL, intL] = label_skel(sk);

pp = regionprops(skL, 'Area','BoundingBox');

for ii = 1:length(pp)
    
    [rr,cc,zz] = getBoundingBox3(pp,ii,3,sz);
    skk = skL(rr,cc,zz)==ii;
    intLL = intL(rr,cc,zz);
    int_neighbors = unique(nonzeros(find_all_neighbors(double(skk), double(intLL), 26)));
    
    sd(ii).intersections = int_neighbors;
    sd(ii).segNeighbors = [];
    sd(ii).label = ii;
    sd(ii).sk = skk;
    sd(ii).intersectionsIm = ~~intLL;
    sd(ii).rr = rr; sd(ii).cc = cc; sd(ii).zz = zz;
    sd(ii).length = pp(ii).Area;
    
    % all segs start out on
    sd(ii).flag = 1;
    
    % no segs are forced on
    sd(ii).hardFlag = 0;
    
end



%----------------------------------------------------------
%
% function to delete spurs of given length
%
%----------------------------------------------------------

function sk = deleteSpur(sk,len)


[sd, skL, intL] = make_skel_data(sk);

if length(sd)==1
    return;
end

for ii = 1:length(sd)
    if sd(ii).length == len
        if length(sd(ii).intersections)<2
            sk(skL==ii) = 0;
        end
    end
end

sk = skel3(double(sk));





