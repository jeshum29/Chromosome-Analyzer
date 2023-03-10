function nd = optiAxis3(nd, OPTZ_METHOD, MAX_AXIS_LENGTH)


if ~exist('OPTZ_METHOD', 'var')
    OPTZ_METHOD = 'intensity';
end

SPUR_LEN = 4;

nd.optParams.OPTZ_METHOD = OPTZ_METHOD;
nd.optParams.MAX_AXIS_LENGTH = MAX_AXIS_LENGTH;

sd = nd.sd;
surfL = nd.surfL;

numS = length(sd);

for ii = 1:numS
    eig_array_z(ii) = nd.sd{ii}.surf_eig_vecs(3,1);
end

eig_array_z = abs(eig_array_z);

for ii = 1:numS
    intensities(ii) = nd.sd{ii}.surface_mean_intensity;
end

intensity_scores = 1 + erf(  (intensities - mean(intensities)) / sqrt(2*(var(intensities))) );
intensity_scores = intensity_scores / max(intensity_scores);

reg = nd.reg;
sz = size(nd.im);

nd.surf_status = zeros(1,numS);
nd.surf_status_init = zeros(1,numS);
nd.surf_status_possible = zeros(1,numS);

%----------------------------------------------------------
%
% find a list of surfaces that separate regions that each
% border more than two regions
%
%----------------------------------------------------------

surfs3 = [];
[regRegArray, surfRegArray, regSurfArray] = makeRegRegArray([], nd.sd, 1:length(nd.ppR), (1:length(nd.ppS)));

for ii = 1:numS
    regs = surfRegArray(ii,:);
    
    if sum(~~regRegArray(regs(1),:)) > 2 && sum(~~regRegArray(regs(2),:)) > 2
        surfs3 = [surfs3, ii];
    end
    
end

%----------------------------------------------------------
%
% SURFACES THAT MUST BE ON
%
%----------------------------------------------------------

% MIN_SURF_RAD = 3;
%
% for ii = 1:numS
%
%     on_flag = ...
%         length(sd{ii}.surfaceless_props) > 1 || ...
%         nd.ppS(ii).Area < MIN_SURF_RAD^2 || ...
%         nd.ppS(ii).Area < (MIN_SURF_RAD+1)^2 && ~isempty(intersect(ii,surfs3)) || ...
%         eig_array_z(ii) > 3/4 && nd.ppS(ii).Area < (MIN_SURF_RAD+1)^2;
%
%     nd.surf_status(ii) = on_flag;
%     nd.surf_status_init(ii) = on_flag;
%     nd.surf_status_possible(ii) = 0;
%
% end


%----------------------------------------------------------
%
% SURFACES THAT CAN POTENTIALLY BE ON
%
%----------------------------------------------------------

for ii = 1:numS
    
    % if surface has to be on, skip
    if nd.surf_status(ii)
        continue;
    end
    
    sk = sd{ii}.sk;
    surf = sd{ii}.mkL==1/2;
    
    cen = sd{ii}.mkRot_cen(sd{ii}.surf_sum>1/2);
    
    if isempty(cen)
        cen = 0;
    end
    
    on_flag = ...
        ... % surfaces with a strong z component
        eig_array_z(ii) > 1/2 || ...
        ... % dim surfaces
        (intensity_scores(ii) < 1/2 && intensities(ii) < mean(intensities)) || ...
        ... % surfaces that don't cut the skeleton
        isempty(nonzeros(surf(~~sk))) || ...
        ... % surfaces that divide regions bordering more than two regions
        ~isempty(intersect(ii, surfs3)) || ...
        ... % anisotropy of merged regions
        mean(cen) > 999 || ...
        ... % small surfaces
        nd.ppS(ii).Area < 15;
    
    
    nd.surf_status_init(ii) = on_flag;
    nd.surf_status_possible(ii) = on_flag;
    
    
end


%----------------------------------------------------------
%
% turn off the surfaces we just turned on if we can do so
% without generating an intersection in the skeleton
%
% try turning them on from most to least intense
%
%----------------------------------------------------------

[surfOn, surfOff] = make_axis_mask(surfL, nd.surf_status_init);
[surfOnPossible, surfOffPossible] = make_axis_mask(surfL, nd.surf_status_possible);

nd.surfOn = surfOn;
nd.surfOff = surfOff;
nd.surfOnPossible = surfOnPossible;
nd.surfOffPossible = surfOffPossible;


sk = makeSkelFromAxisMask(nd.reg + surfOff, SPUR_LEN);
int = find_intersections_3(double(sk), 26);
num_int = sum(int(:)>2);
turn_surf_off = 0*nd.surf_status;


% list of possible surfaces excluding definite surfaces
inds = find(nd.surf_status_possible .* (~nd.surf_status));

if strcmpi(OPTZ_METHOD, 'intensity')
    
    order_param = intensities(inds);
    [order_param, ord] = sort(order_param,'descend');
    
elseif strcmpi(OPTZ_METHOD, 'orientation')
    
    order_param = eig_array_z(inds);
    [order_param, ord] = sort(order_param,'ascend');
    
elseif strcmpi(OPTZ_METHOD, 'random')
    
    order_param = (1 - eig_array_z) + intensity_scores;
    order_param = order_param(inds);
    
    [order_param, ord] = sort(order_param, 'descend');
    
end

inds = inds(ord);
sz = size(surfL);

pp_sk = regionprops(sk,'Area');
areas_init = [pp_sk.Area];

skels_too_long = find(areas_init > MAX_AXIS_LENGTH);

 
for ii = 1:length(inds)
    
    ind = inds(ii);
    
    [rr,cc,zz] = getBoundingBox3(nd.ppS, ind, 1, sz);
    
    surff = surfL(rr,cc,zz)==ind;
    surfOffTemp = surfOff;
    
    surfOffTemp_ = surfOffTemp(rr,cc,zz);
    surfOffTemp_(surff) = 1;
    surfOffTemp(rr,cc,zz) = surfOffTemp_;
    
    sk = makeSkelFromAxisMask(nd.reg + surfOffTemp, SPUR_LEN);
    
    %     if ind==17 || ind==2
    %        vizVolume2(sk,0,[1,1,1],66,1,1);
    %         '';
    %     end
    
    pp_sk = regionprops(sk,'Area');
    areas = [pp_sk.Area];
    
    areas = setxor(areas, areas_init(skels_too_long));
    int = find_intersections_3(double(sk), 26);
    num_int_new = sum((int(:)>2));
    
    if num_int  == num_int_new && max(areas) < MAX_AXIS_LENGTH
        
        surfOn_ = surfOn(rr,cc,zz);
        surfOff_ = surfOff(rr,cc,zz);
        
        surfOn_(surff) = 0;
        surfOff_(surff) = 1;
        
        surfOn(rr,cc,zz) = surfOn_;
        surfOff(rr,cc,zz) = surfOff_;
        
        nd.surf_status(ind) = 0;
        inds(ii) = 0;
        
    else
        nd.surf_status(ind) = 1;
    end
    
end


[surfOn, surfOff] = make_axis_mask(surfL, nd.surf_status);
[surfOnPossible, surfOffPossible] = make_axis_mask(surfL, nd.surf_status_possible);

mk = nd.reg + surfOff;
mk = bwlabeln(mk, 6);
num_reg = max(mk(:));

% keep turning off segments if there are more than six axes

if num_reg > 6
    
    inds = inds(~~inds);
    
    for ii = 1:length(inds)
        
        ind = inds(ii);
        
%         if ~(nd.surf_status_init(ind)==1 && nd.surf_status(ind)==0)
%             continue;
%         end
        
        [rr,cc,zz] = getBoundingBox3(nd.ppS, ind, 1, sz);
        
        surff = surfL(rr,cc,zz)==ind;
        surfOffTemp = surfOff;
        
        surfOffTemp_ = surfOffTemp(rr,cc,zz);
        surfOffTemp_(surff) = 1;
        surfOffTemp(rr,cc,zz) = surfOffTemp_;
        
        mk = nd.reg + surfOffTemp;
        mk = bwlabeln(mk, 6);
        num_reg_new = max(mk(:));
        
        if (num_reg_new < num_reg)
            
            surfOn_ = surfOn(rr,cc,zz);
            surfOff_ = surfOff(rr,cc,zz);
            
            surfOn_(surff) = 0;
            surfOff_(surff) = 1;
            
            surfOn(rr,cc,zz) = surfOn_;
            surfOff(rr,cc,zz) = surfOff_;
            
            nd.surf_status(ind) = 0;
            num_reg = num_reg_new;

        else
            nd.surf_status(ind) = 1;
        end
        
        if num_reg==6
            break;
        end
        
    end
end


[surfOn, surfOff] = make_axis_mask(surfL, nd.surf_status);
[surfOnPossible, surfOffPossible] = make_axis_mask(surfL, nd.surf_status_possible);

nd.surfOn = surfOn;
nd.surfOff = surfOff;
nd.surfOnPossible = surfOnPossible;
nd.surfOffPossible = surfOffPossible;

nd.ordered_surf_inds_for_optz = inds;

nd.surf_intensities = intensities;
nd.surf_eig_array_z = eig_array_z;



return;

'';


















%----------------------------------------------------------
%
% make region-region and region-surface arrays
%
%----------------------------------------------------------


[regRegArrayRaw,...
    surfRegArrayRaw,...
    regSurfArrayRaw] = makeRegRegArray([], nd.sd, 1:length(nd.ppR), 1:length(nd.ppS));

groupsRaw = makeSubGroups(regRegArrayRaw, 1:length(nd.ppR));


[regRegArray,...
    surfRegArray,...
    regSurfArray] = makeRegRegArray([], nd.sd, 1:length(nd.ppR), setxor(1:length(nd.ppS), find(nd.surf_status_init)));

groups = makeSubGroups(regRegArray, 1:length(nd.ppR));


%----------------------------------------------------------
%
% another way to make these arrays involves generating
% a mask and bwlabeling. these masks are useful for other
% reasons
%
%----------------------------------------------------------

mkRaw = bwlabeln(nd.reg + ~~nd.surfL, 6);
ppRaw = regionprops(mkRaw, 'Area');

mk = bwlabeln(nd.reg + nd.surfOff, 6);
pp = regionprops(mk, 'Area');

%----------------------------------------------------------
%
% make a skeleton with the same labeling as the mask
%
%----------------------------------------------------------

sk = 0*mk;

for ii = 1:length(pp)
    
    skk = makeSkelFromAxisMask(double(mk==ii), SPUR_LEN);
    sk = sk + (~~skk)*ii;
    
end

ppSK = regionprops(sk,'Area');


%----------------------------------------------------------
%
% hard thresh to find axes that are the right length
%
%----------------------------------------------------------

areas = [ppSK.Area];
probs = calcAxisProb([ppSK.Area]);
REAL_INDS = find( (probs > PROB_THRESH) .* (probs < 1 - PROB_THRESH) );

%----------------------------------------------------------
%
% include unusually short axes that cannot be made any longer
% but turning off surfaces
%
%----------------------------------------------------------

for ii = 1:length(ppSK)
    
    if probs(ii) < PROB_THRESH && areas(ii) > 30
        
        % search for a raw group that exactly matches
        % this axis's subgroup. if one exists, the axis
        % cannot be made any longer
        
        for jj = 1:length(ppRaw)
            
            if isempty(setdiff(groups(ii).regs, groupsRaw(jj).regs))
                
                REAL_INDS = [REAL_INDS, ii];
                
            end
        end
    end
end

NUM_REAL = length(REAL_INDS);

%----------------------------------------------------------
%
% finally, the constraint that there must be 6 axes
%
%----------------------------------------------------------

NUM_MISSING_AXES = 6 - NUM_REAL;

LONG_INDS = find(probs > 1 - PROB_THRESH);

% vizNucleusSIM_config(mk==LONG_INDS(1), nd.surfL, regs, surfs_off, [surfs_on, surfs(configs(jj,:))], nd.ppS,[]);


% if there's just one overseg'd region, it's easier
if length(LONG_INDS)==1
    
    regs = groups(LONG_INDS(1)).regs;
    
    % find surfs to consider turning on to break the region up
    surfs_all = nonzeros(unique(regSurfArrayRaw(regs, :)));
    
    surfs = intersect(surfs_all, surfs3);
    
    surfs_on = intersect(surfs_all, nd.good_surfs);
    
    % don't consider surfaces that are already on
    surfs = setdiff(surfs, surfs_on);
    
    N = length(surfs);
    numPerms = 2^N
    
    for num_on = 1:N
        
        num_on
        configs = logical(unique(perms([ones(1, num_on), zeros(1, N - num_on)]),'rows'));
        
        for jj = 1:size(configs,1)
            
            surfs_off = setdiff(surfs_all, [surfs_on, surfs(configs(jj, :))]);
            
            [rrA, srA, rsA] = makeRegRegArray([], nd.sd, regs, surfs_off);
            
            subGroups = makeSubGroups(rrA, regs);
            
            
            if length(subGroups) == NUM_MISSING_AXES
                
                sum([nd.ppR(subGroups(1).regs).Area])
                sum([nd.ppR(subGroups(2).regs).Area])
                
                vizNucleusSIM_config(mk==LONG_INDS(1), nd.surfL, regs, surfs_off, [surfs_on, surfs(configs(jj,:))], nd.ppS,[]);
                
                '';
                
            end
            
        end
    end
end




'';


end

function [surfOn, surfOff] = make_axis_mask(surfL, surf_status)

surfOn = 0*surfL;
surfOff = 0*surfL;

for ii = 1:length(surf_status)
    
    if surf_status(ii)==10
        surfOn = surfOn + (surfL==ii)*2;
    elseif surf_status(ii)==1
        surfOn = surfOn + (surfL==ii);
    else
        surfOff = surfOff + (surfL==ii);
    end
end

end

