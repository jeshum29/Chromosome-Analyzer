function nd = optiAxis3Smart(nd, dispFlag, isOMX)


% ============================================================
%
% Function to iteratively optimization axis segmention by choosing which
% surfaces should segment the image and which are spurious
%
% Keith Cheveralls
% Dernburg Lab
% March 2013
%
% May 2014: Updated to accomodate OMX images for which the training
% data set is relatively sparse and the calcSurfGoodProb function
% is not used
%
% ============================================================
    

if ~exist('isOMX', 'var')
    isOMX = 0;
end

% ============================================================
%
%  Hard-coded parameters
%
% ============================================================


SPUR_LEN               = 3;
SKEL_CROP_RADIUS       = 20;
ABS_MIN_AXIS_VOL       = 3000;
BAD_SURF_MAX_PROB      = 0.02;
GOOD_SURF_MIN_PROB     = 0.98;
MIN_SURF_AREA          = 5;

MAX_AXIS_LENGTH_OMX = 160;

if ~exist('dispFlag')
    dispFlag = 0;
end

if isempty(nd.sd)
    return;
end

% ============================================================
%
%  Initialization stuff
%
% ============================================================

nd.optiData.status                = [];
nd.optiData.axis_vol_score        = [];
nd.optiData.longest_axis_score    = [];
nd.optiData.problem               = [];
nd.optiData.status_explicit       = [];
nd.optiData.success_flag          = [];

nd.focusGroups = [];

sz        = size(nd.surfL);
surfL     = nd.surfL;
reg       = nd.reg;
regProps  = regionprops(nd.regL, 'Centroid');
numReg    = length(nd.ppR);

FOCI_PRESENT = 0;
if isfield(nd, 'fd')
    if ~isempty(nd.fd)
        FOCI_PRESENT = 1;
    end
end

if isOMX
    FOCI_PRESENT = 0;
end


% ============================================================
%
%  Calculation of surface scores and probabilities
%
% ============================================================

if ~isOMX
    % load calculated score parameters and probability parameters
    [scoreParams, probGoodParams] = loadScoreFunctionParams3;

    % load the data set used to calculate those parameters
    dataGlobal = loadGlobalData3;
else
    scoreParams = loadScoreFunctionParamsOMX;
    dataGlobal = loadGlobalDataOMX;
end

% make surface properties for this nucleus
[props, status] = makePropsArray2v3(nd);
data(1).props = props;
data(1).status = status;

% calculate surface properties for this nucleus (requires the training data)
avProps = calcAvSurfProps3(data, dataGlobal);

% now calculate scores
scores = calcSurfScoreFunction3(avProps, scoreParams, isOMX);

if ~isOMX
    % and the probabilities that each surface is good
    probGood = calcSurfGoodProb(scores, probGoodParams);
else
    probGood = 0*scores + .5;
end
    

% ============================================================
%
%  Calculation of various axis-length distributions from the training data set
%
% ============================================================

if ~isOMX
    % here we calculate the mean and inverse covariance
    % for the set of five normalized axis lengths
    % usage: % exp( - 1/2 * (len - gauss.M) * gauss.C * (len - gauss.M)');
    gauss_axis_len = calcAxisLengthDistro3(dataGlobal);

    % finally we compute the same quantities for the shortest and longest axis
    [longestAxisDistro, shortestAxisDistro] = calcAxisLengthDistroRel3(dataGlobal);
else
    
end

% ============================================================
%
%  Establish initial surface statuses
%
% status = 1 means the surface DOES segment the image
%
% ============================================================

% start with all surfs on (i.e., all regions separated)
status = 0*scores + 1;

% merge regions separated by surfaces that are almost definitely spurious
status(probGood < BAD_SURF_MAX_PROB) = 0;

% force very small surfaces and high-probability surfaces on
status_force_on = 0*scores;
status_force_on([nd.ppS.Area] < MIN_SURF_AREA) = 1;
status_force_on(probGood > GOOD_SURF_MIN_PROB) = 1;

disp = 0;

if disp
    figNum = 1;
    [surfOn, surfOff] = make_axis_mask(nd.surfL, status);
    mk = reg + surfOff;
    mk = ~~mk;
    [mkL, num] = bwlabeln(mk, 6);
    cc = colormap(lines(max(mkL(:))));
    vizSegColor(mkL, figNum, 1, cc);
end



% ============================================================
%
%  OPTIMIZATION PREP PART 1
%
%  Make region groups based on the set of surfaces that are hard-coded on
%
% ============================================================

[surfOn, surfOff] = make_axis_mask(surfL, status);

regRegArray = makeRegRegArray([], nd.sd, 1:numReg, find(~status));
groups = makeSubGroups(regRegArray, 1:numReg);

% ============================================================
%
%  OPTIMIZATION PREP PART 2
%
%  Make an initial skeleton and  calculate the total axis length
%  This will be used later to calculate the probability that an axis is too long
%
% ============================================================

sk_max = make_skel(~~(reg + surfOff + surfOn), 0);

for ii = 1:length(nd.ppR)
    len = sum(sk_max(nd.regL(:)==ii));
    nd.ppR(ii).skelLength = len;
end

for jj = 1:length(groups)
    axisLengths(jj) = sum([nd.ppR(groups(jj).regs).skelLength]);
end

max_len = max(axisLengths);
totalAxisLength = sum([nd.ppR.skelLength]);


% ============================================================
%
%  OPTIMIZATION PREP PART 3
%
%  If foci are present, determine the region nearest to each focus
%  For now, we assume that each focus labels a different chromosome
%
% ============================================================

if FOCI_PRESENT
    
    fd = nd.fd;
    
    % find the region nearest to each focus
    for nn = 1:length(fd)
        
        pos = fd{nn}.position;
        
        for ii = 1:numReg
            dd(ii) = norm(regProps(ii).Centroid - pos([2,1,3]));
        end
        
        [dd_, ord] = sort(dd, 'ascend');
        reg_id = ord(1);
               
        % the following code may cause more errors than it prevents
        
        % check to see if a focus has already been associated with this region
%         if nn>1
%             if ~isempty(intersect(focusRegions, ord(1)))
%                 if isempty(intersect(focusRegions, ord(2)))
%                     reg_id = ord(2);
%                 end
%             end
%         end
        
        focusRegions(nn) = reg_id;
        fd{nn}.region = reg_id;
    end
    
    % determine if any foci share a region group - this is bad news since
    % it means some of the surfaces we've already turned off should be on
    sharedFocus = 0*focusRegions;
    
    for nn = 1:length(fd)
        for jj = 1:length(groups)
            if ~isempty(intersect(groups(jj).regs, fd{nn}.region))
                fd{nn}.groupID_current = jj;
                break;
            end
        end
    end
    
    for nn = 1:length(fd)
        for nnn = nn+1:length(fd)
            if fd{nn}.groupID_current==fd{nnn}.groupID_current
                sharedFocus(nnn) = nn;
                break;
            end
        end
    end
    
    for nn = 1:length(fd)
        fd{nn}.sharedFocusIndex = sharedFocus(nn);
    end
    
end


% find the set of merged regions that coincides with each focus

if FOCI_PRESENT
    for nn = 1:length(fd)
        
        if fd{nn}.sharedFocusIndex
            fd{nn}.regs = [];
            continue;
        end
        
        for ii = 1:length(groups)
            if ~isempty(intersect(groups(ii).regs, fd{nn}.region))
                fd{nn}.regs = groups(ii).regs;
                break;
            end
        end
        
    end
    
    
    nd.fd = fd;
    
end

% ============================================================
%
%  Iterative optimization: STEP #1
%
% Look at each surface in probability order.
%
% Turn it *off* if doing so does not 
% 1) generate an intersection
% 2) generate a too-long axis
% 3) merge axes corresponding to different FISH foci
%
% Stop when there are 6 regions.
%
% ============================================================

inds = find(status.*(~status_force_on));
[scores_,ord] = sort(scores(inds), 'ascend');

inds = inds(ord);
tic;

for ii = 1:length(inds)
    
    % ============================================================
    %
    %  Crop around the surface and skeletonize the cropped region to
    %  determine whether the surface will generate a new intersection if
    %  turned off
    %
    % ============================================================
    
    [rr,cc,zz] =    getBoundingBox3(nd.ppS, inds(ii), 1, sz);
    [rrW,ccW,zzW] = getBoundingBox3(nd.ppS, inds(ii), SKEL_CROP_RADIUS, sz);
    
    surff = surfL(rr,cc,zz)==inds(ii);
    
    surfOffTemp =               surfOff;
    surfOffTemp_crop =          surfOffTemp(rr,cc,zz);
    surfOffTemp_crop(surff) =   1;
    surfOffTemp(rr,cc,zz) =     surfOffTemp_crop;
    
    sk_old = make_skel( reg(rrW,ccW,zzW) + surfOff(rrW,ccW,zzW),     SPUR_LEN);
    sk_new = make_skel( reg(rrW,ccW,zzW) + surfOffTemp(rrW,ccW,zzW), SPUR_LEN);
    
    num_int_old = count_intersections(sk_old);
    num_int_new = count_intersections(sk_new);
    
    % ============================================================
    %
    %  Use length distributions to make sure the surface won't generate an
    %  unreasonably long axis
    %
    % ============================================================

    status(inds(ii)) = 0;
    regRegArray = makeRegRegArray([], nd.sd, 1:numReg, find(~status));
    groups = makeSubGroups(regRegArray, 1:numReg);
    
    clear axisLengths
    for jj = 1:length(groups)
        axisLengths(jj) = sum([nd.ppR(groups(jj).regs).skelLength]);
    end
    
    if ~isOMX
        vec = [max(axisLengths), totalAxisLength] - longestAxisDistro.M;
        score_max = (vec * longestAxisDistro.C * vec') * (max(axisLengths) > longestAxisDistro.M(1));
    end
    
    % ============================================================
    %
    %  If foci are present, determine if turning off the region would merge the regions
    % associated with two different FISH foci
    %
    % ============================================================
    
    focusConflictFlag = 0;
    
    if FOCI_PRESENT
        
        for nn = 1:length(fd)
            
            if fd{nn}.sharedFocusIndex
                fd{nn}.groupID_current = 0;
                continue;
            end
            
            for jj = 1:length(groups)
                if ~isempty(intersect(groups(jj).regs, fd{nn}.region))
                    fd{nn}.groupID_current = jj;
                    break;
                end
            end
            
        end
        
        for nn = 1:length(fd)
            
            if ~fd{nn}.groupID_current
                continue;
            end
            
            for nnn = nn+1:length(fd)
                if fd{nn}.groupID_current==fd{nnn}.groupID_current
                    focusConflictFlag = 1;
                end
            end
            
        end
        
    end
    
    % ============================================================
    %
    %  Here is the decision about whether to turn off the surface
    %
    if ~isOMX
        
        turn_off = ...
            (num_int_new <= num_int_old) &&...
            (score_max < 10 || max(axisLengths)==max_len) &&...
            focusConflictFlag==0;
        
    else
        
        turn_off = ...
            (num_int_new <= num_int_old) &&...
            (max(axisLengths) < MAX_AXIS_LENGTH_OMX) &&...
            focusConflictFlag==0;
        
    end
    %
    % ============================================================
    
    if 0
        
        figNum = 11;
        [surfOn_display, surfOff_display] = make_axis_mask(nd.surfL, status);
        mk = ~~(reg + surfOff_display);
        [mkL, num] = bwlabeln(mk, 6);
        cmap = colormap(lines(max(mkL(:))));
        vizSegColor(mkL, figNum, 1, cmap);
        
        figNum = 22;
        vizMask(~~surfOn,           [1,0,0],   figNum, 1, 1);
        vizMask(~~surfOff,          [0,1,0],   figNum, 0, 0);
        vizVolume3(surfL==inds(ii), [1,1,0],   figNum, 0, 0);
        vizMask(~~nd.regL,          [1,1,1]/2, figNum, 0, 0, 1);
        
        '';
        
    end
    
    
    if turn_off
        
        surfOn_crop = surfOn(rr,cc,zz);
        surfOff_crop = surfOff(rr,cc,zz);
        
        surfOn_crop(surff) = 0;
        surfOff_crop(surff) = 1;
        
        surfOn(rr,cc,zz) = surfOn_crop;
        surfOff(rr,cc,zz) = surfOff_crop;
        
        max_len = max(axisLengths);
        
    else
        
        status(inds(ii)) = 1;
        
    end
    
    
end

toc;

% [surfOn, surfOff] = make_axis_mask(nd.surfL, status);
% mk = ~~(reg + surfOff);
% [mkL, num] = bwlabeln(mk, 6);

nd.optiData.status =                      status;










return;











if dispFlag
    
    figNum = 1;
    cmap = colormap(lines(max(mkL(:))));
    vizSegColor(mkL, figNum, 9, cmap);
    
end


pp = regionprops(mkL,'Area');
a = [pp.Area];
a_raw = a;
vec = [max(a), nukeVolume] - longestAxisDistro.M;
longest_axis_score = (vec * longestAxisDistro.C * vec') * (max(a) > longestAxisDistro.M(1));


if length(a)>5
    a = a/sum(a);
    a = sort(a,'descend');
    a = a(1:5);
    
    axis_vol_score = ( (a - gauss_axis_vols.M) * gauss_axis_vols.C * (a - gauss_axis_vols.M)' );
else
    axis_vol_score = NaN;
end

AXIS_VOL_SCORE_MAX = 30;
LONGEST_AXIS_SCORE_MAX = 10;

underSegmentation = -1;
overSegmentation = 1;
noProblem = 0;

problem = noProblem;

% if the axis score is too high and the longest axis was too long
% this is an underseg problem -- we turned too many surfs off
if axis_vol_score > AXIS_VOL_SCORE_MAX && longest_axis_score > LONGEST_AXIS_SCORE_MAX
    
    problem = underSegmentation;
    
    % if the axis score is too high but the longest axis is not too long...
    % this suggests oversegmentation -- we didn't turn enough surfs off
elseif axis_vol_score > AXIS_VOL_SCORE_MAX && longest_axis_score <= LONGEST_AXIS_SCORE_MAX
    
    problem = overSegmentation;
    
    % if there were less than six regions...
    % this is definitely an underseg problem
elseif isnan(axis_vol_score)
    
    problem = underSegmentation;
    
    % if the axis score was okay but one axis is still too long...
    % this may not be a problem
elseif axis_vol_score <= AXIS_VOL_SCORE_MAX && longest_axis_score > LONGEST_AXIS_SCORE_MAX
    
    problem = noProblem;
    
end

status_explicit = [];
success_flag = 0;

switch problem
    
    case underSegmentation
        
        try
            [status_explicit, success_flag] = optiAxis3SmartExplicit2v2(nd);
        catch
            'Explicit optimization error'
            nd.dir
        end
        
    case overSegmentation
        
        try
            [status_explicit, success_flag] = optiAxis3SmartExplicit2v2(nd);
        catch
            'Explicit optimization error'
            nd.dir
        end
        
    case noProblem
        
        status_explicit = [];
        success_flag = [];
end

nd.optiData.status =                 status;
nd.optiData.axis_vol_score =         axis_vol_score;
nd.optiData.longest_axis_score =     longest_axis_score;
nd.optiData.problem =                problem;
nd.optiData.status_explicit =        status_explicit;
nd.optiData.success_flag =           success_flag;




'';

% compute axis volumes for a given mask volume so that we can choose
% the set of good surfs that generates axis lengths closest to the
% average length for the nucleus size

% this will avoid having to do any skeletonization at all

% also, if we reduce to set of region connectivities, we can avoid any
% image manipulation entirely

% write code to run the original opti3AxisQuick without altering any other
% fields in the nd files

end

function [surfOn, surfOff] = make_axis_mask(surfL, surf_status)

surfOn = 0*surfL;
surfOff = 0*surfL;

for ii = 1:length(surf_status)
    
    if surf_status(ii)==1
        surfOn(surfL==ii) = 1;
    else
        surfOff(surfL==ii) = 1;
    end
    
end

end

function sk = delete_spurs(sk, num)

for ii = 1:num
    
    spurs = count_neighbor_pixels(double(sk), 26)==1;
    sk(spurs) = 0;
    
end

end


function sk = make_skel(mk, spur_len)

mk = double(imfill(~~mk,'holes'));

sk = skel3(double(~~mk));
sk = delete_spurs(sk, spur_len);

int = count_neighbor_pixels(double(sk), 26)>2;
sk(int) = 0;

% delete isolated pixels
sk = ~~count_neighbor_pixels(double(sk), 26);

sk(int) = 1;

sk = skel3(double(sk));

end

function n = count_intersections(sk)

int = find_intersections_3(double(sk), 26);
n = length(find(int(:)>2));

end

