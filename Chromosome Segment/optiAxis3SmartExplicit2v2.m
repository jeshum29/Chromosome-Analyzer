function [status, success_flag] = optiAxis3Smart(nd, dispFlag)

SPUR_LEN = 3;
SKEL_CROP_RADIUS = 20;
ABS_MIN_AXIS_VOL = 3000;
BAD_SURF_MAX_PROB = .03;
GOOD_SURF_MIN_PROB = .98;
MIN_SURF_AREA = 5;

N_MAX = 14;

success_flag = 1;
status = [];

if ~exist('dispFlag')
    dispFlag = 0;
end

% load calculated score parameters and probability parameters
[scoreParams, probGoodParams] = loadScoreFunctionParams;

% load the data set used to calculate those parameters
% i.e., the training data set
dataGlobal = loadGlobalData;

% make surface properties for this nucleus
props = makePropsArray2v2(nd);
data(1).props = props;

% calculate surface properties for this nucleus using surface data
% and certain averages from the training data set
[props, status] = calcSurfProps(data, dataGlobal);

% now we can calculate scores
scores = calcSurfScoreFunction(props, scoreParams);

% and the probabilities that each surface is good
probGood = calcSurfGoodProb(scores, probGoodParams);


sz = size(nd.surfL);
surfL = nd.surfL;
reg = nd.reg;
numReg = length(nd.ppR);
nuke_vol = data(1).props(9,1);
areas = [nd.ppS.Area];

% here we calculate the mean and inverse covariance
% for the set of five normalized axis volumes
gauss_axis_vols = calcAxisVolumeDistro(dataGlobal);

% exp( - 1/2 * (vol - gauss.M) * gauss.C * (vol - gauss.M)');

% finally we compute the probability for the shortest and longest axis
% given the nucleus volume/radius

[gauss_longest, gauss_shortest] = calcAxisVolumeNukeVolumeDistro(dataGlobal);

a = 0:100:1e5;
for ii = 1:length(a)
    vec = [a(ii), nuke_vol] - gauss_longest.M;
    prob(ii) = exp(-(vec * gauss_longest.C * vec') * (max(a) > gauss_longest.M(1)));
end

av_max_vol = sum(prob.*a)/sum(prob);

% status = 1 means the surface IS present in the image

% start with all surfs on (i.e., all regions separated)

status = 0*scores + 1;

% turn off surfaces that are almost definitely bad
status(probGood < BAD_SURF_MAX_PROB) = 0;

status_force_on = 0*scores;

% these surfaces are forced on
status_force_on(areas < MIN_SURF_AREA) = 1;
status_force_on(probGood > GOOD_SURF_MIN_PROB) = 1;


always_on = status_force_on;
always_off = ~status;


disp = 0;

if disp
    figNum = 1;
    [surfOn, surfOff] = make_axis_mask(nd.surfL, status);
    mk = reg + surfOff;
    mk = ~~mk;
    [mkL, num] = bwlabeln(mk, 6);
    cc = colormap(lines(max(mkL(:))));
    vizSegColor(mkL, figNum, 1, cc);
    
    vizSegInitial_optiAxisExplicit(nd, 7, 1, always_off, always_on);

    
end

% first step: look at each surface in probability order (most likely to be bad, i.e., should be off)
% and turn it *off* if doing so does not generate an intersection or a very large region.
% Stop when there are 6 regions.

inds = find(status.*(~status_force_on));
[scores_,ord] = sort(scores(inds), 'ascend');

inds = inds(ord);

[surfOn, surfOff] = make_axis_mask(surfL, status);

regRegArray = makeRegRegArray([], nd.sd, 1:numReg, find(~status));
groups = makeSubGroups(regRegArray, 1:numReg);

for ii = 1:length(groups)
    areas(ii) = sum([nd.ppR(groups(ii).regs).Area]);
end

regions = 1:numReg;
surfaces = 1:length(status);
surfaceData = nd.sd;


for ii = 1:length(surfaces)
    
    nnbs = nonzeros(unique(surfaceData{ii}.neighbors));
    neighborsArray(ii,1:length(nnbs)) = nnbs;
    
end


tic;

N = length(inds);
clear areas;

if N < N_MAX
    
    'Entering explicit optimization!'
    ['N = ' num2str(N)]
    
    for nn = 0:2^N-1
        
        sub_status = sprintf(['%0' num2str(N) 'ld'], str2num(dec2base(nn, 2)));
        
        for ii = 1:N
            status(inds(ii)) = str2double(sub_status(ii));
        end
        
%         regRegArray = make_regreg_array(regionData, surfaceData, 1:numReg, find(~status));
%         groups_ = makeSubGroups(regRegArray, 1:numReg);
        
        neighborsArray_ = neighborsArray;
        neighborsArray_(find(status),:) = [];
        groups = make_groups(neighborsArray_,1:numReg);
        

                for ii = 1:size(groups,1)
                    areas(nn+1,ii) = sum([nd.ppR(nonzeros(groups(ii,:))).Area]);
                    log_prob(nn+1,ii) = sum(log(probGood(~~status)));
                    num_surf_on(nn+1,ii) = sum(~~status);
                end

    end
    
else
    
    success_flag = 0;
    return;
    
end

toc;

vscores = batchCalcAxisVolScores(areas);
[v,ord] = sort(vscores);

dataGlobal = loadGlobalData;
gauss_axis_vols = calcAxisVolumeDistro(dataGlobal);
[gauss_longest, gauss_shortest] = calcAxisVolumeNukeVolumeDistro(dataGlobal);

nuke_vol = sum([nd.ppR.Area]);

for ii = 1:size(areas,1)
    vec = [max(a), nuke_vol] - gauss_longest.M;
    mscores(ii) = (vec * gauss_longest.C * vec') * (max(a) > gauss_longest.M(1));
end

nReg = sum(areas>ABS_MIN_AXIS_VOL/5,2);
possible_configs = (nReg==6);

total_score = mscores/max(mscores) + vscores/max(vscores);

total_score(~possible_configs) = 1e9;

[t,ord] = sort(total_score);

sub_status = sprintf(['%0' num2str(N) 'ld'], str2num(dec2base(ord(1)-1, 2)));
        
for ii = 1:N
    status(inds(ii)) = str2double(sub_status(ii));
end


if dispFlag
    
    figNum = 11;
    [surfOn, surfOff] = make_axis_mask(nd.surfL, status);
    mk = ~~(reg + surfOff);
    [mkL, num] = bwlabeln(mk, 6);
    cmap = colormap(lines(max(mkL(:))));
    vizSegColor(mkL, figNum, 9, cmap);
    
end


'';

end

function groups = make_groups(neighborsArray, regs_all)

sz = size(neighborsArray);

regs = nonzeros(unique(neighborsArray(:)));
isolatedRegs = (setxor(regs, regs_all));

if ~isempty(isolatedRegs)
    cc = length(isolatedRegs);
    groups(1:cc) = isolatedRegs;
    groups = groups';
else
    cc = 0;
end


for ii = 1:length(regs)
    
    [rows,cols] = find(neighborsArray==regs(ii));
    
    if isempty(rows)
        continue;
    end
    
    cc = cc+1;
    group_regs = [];
    go = 1;
    
    while go
        regs_ = nonzeros(unique(neighborsArray(rows,:)))';
        
        if ~isempty(regs_)
            neighborsArray(rows,:) = [];
            rows = [];
            
            for jj = 1:length(regs_)
                [rows_, cols_] = find(neighborsArray==regs_(jj));
                
                if ~isempty(rows_)
                    rows = [rows, rows_'];
                end
                
            end
            
            group_regs = unique([group_regs, regs_]);
            
        else
            rows = [];
        end
        
        if isempty(rows)
            go = 0;
        end
        
    end
    
    groups(cc,1:length(group_regs)) = group_regs;
    
end


end





%
%
% [surfOn, surfOff] = make_axis_mask(nd.surfL, status);
% mk = ~~(reg + surfOff);
% [mkL, num] = bwlabeln(mk, 6);
% pp = regionprops(mkL,'Area');
% a = [pp.Area];
% a_raw = a;
% vec = [max(a), nuke_vol] - gauss_longest.M;
% longest_axis_score = (vec * gauss_longest.C * vec') * (max(a) > gauss_longest.M(1));
%
%
% if length(a)>5
%
%     a = a/sum(a);
%     a = sort(a,'descend');
%     a = a(1:5);
%
%     axis_vol_score = ( (a - gauss_axis_vols.M) * gauss_axis_vols.C * (a - gauss_axis_vols.M)' );
% else
%     axis_vol_score = NaN;
% end
%
% AXIS_VOL_SCORE_MAX = 30;
% LONGEST_AXIS_SCORE_MAX = 10;
%
% underSegmentation = -1;
% overSegmentation = 1;
% goodSegmentation = 0;
%
% problem = goodSegmentation;
%
%
%
%
%
% nd.surfOn = surfOn;
% nd.surfOff = surfOff;
% nd.surf_status = status;
% nd.axis_vol_score = axis_vol_score;
%
%
%
% disp = 0;
% if disp
%     figNum = 1;
%     cmap = colormap(lines(max(mkL(:))));
%     vizSegColor(mkL, figNum, 9, cmap);
%
% end
%
% '';



% compute axis volumes for a given mask volume so that we can choose
% the set of good surfs that generates axis lengths closest to the
% average length for the nucleus size

% this will avoid having to do any skeletonization at all

% also, if we reduce to set of region connectivities, we can avoid any
% image manipulation entirely



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

