function [mkF, props] = makeNukesFromMask(bandpass_image, GONAD)

%-------------------------------------------------
%
% OPERATES ON FULL SIZED MASK
%
%-------------------------------------------------

MIN_AXIS_VOLUME = GONAD.MIN_AXIS_VOLUME;
MAX_AXIS_VOLUME = GONAD.MAX_AXIS_VOLUME;
MAX_NUKE_DIA = GONAD.MAX_NUKE_DIA;
MIN_NUKE_DIA = GONAD.MIN_NUKE_DIA;

ASPECT_RATIO = GONAD.ASPECT_RATIO;

mkL = bwlabeln(bandpass_image > GONAD.LOW_THRESH, 6);
mkL_high = bandpass_image > GONAD.HIGH_THRESH;

props = regionprops(mkL, 'Area','Centroid','BoundingBox');
sz = size(mkL);

for ii = 1:length(props)
    
    if props(ii).Area < MIN_AXIS_VOLUME
        
        [r,c,z] = getBoundingBox3(props,ii,3,sz);
        mkk = mkL(r,c,z);
        mkk(mkk==ii) = 0;
        mkL(r,c,z) = mkk;
        
    end
end

mkL = bwlabeln(~~mkL, 6);
props = regionprops(mkL, 'Area','Centroid','BoundingBox');

for ii = 1:length(props)
    
    include_flag = isNuke(props, ii, 0, MAX_AXIS_VOLUME, MAX_NUKE_DIA, sz);

    if ~include_flag
        
        [r,c,z] = getBoundingBox3(props,ii,3,sz);
        
        mkk = mkL(r,c,z)==ii;
        mkk_high = mkL_high(r,c,z);
        mkk = mkk.*mkk_high;
        
        mkk_copy = mkL(r,c,z);
        mkk_copy(mkL(r,c,z)==ii) = mkk(mkL(r,c,z)==ii)*ii;
        
        mkL(r,c,z) = mkk_copy;
        
    end
    
end

mkL = bwlabeln(~~mkL, 6);
props = regionprops(mkL, 'Area','Centroid','BoundingBox');

numNukes = length(props);
dists = [20:20:MAX_NUKE_DIA+40];
mkF = mkL;

for ii = 1:length(dists)
    
    mkF = merge_regions(mkF, props, dists(ii), MIN_AXIS_VOLUME, MAX_AXIS_VOLUME, MAX_NUKE_DIA, ASPECT_RATIO, sz);
    props = regionprops(mkF, 'Area','Centroid','BoundingBox');
    
    
    numNukes = length(props);
    
end

%  save3(uint16(mkF), [GONAD.writeDir filesep], 'mkF3', 0);

% dd = [gonad.writeDir slash];

% save3(uint16(mkL), dd, 'mkBL', 0);
% save3(uint16(mkF), dd, 'mkBF', 0);
% save3(autogain(mkF), [dd slash 'mkF'], 'mkF', 1);
% save3c(colorMask3((mkF)), dd, 'mkF_color', 0);
% gonad.props = props;
% save([dd slash 'gonad.mat'], 'gonad');

end

function mkF = merge_regions(mkL, props, DIST_THRESH, MIN_AXIS_VOLUME, MAX_AXIS_VOLUME, MAX_NUKE_DIA, ASPECT_RATIO, sz)

for ii = 1:length(props)
    props(ii).neighbors = [];
    props(ii).dists = [];
    
    for jj = 1:length(props)
        
        if ii==jj
            continue;
        end
        
        dd = props(ii).Centroid - props(jj).Centroid;
        
        dd = [dd(1), dd(2), dd(3)*ASPECT_RATIO];
        
        dd = sqrt(sum(dd.^2));
        
        if dd < DIST_THRESH
            props(ii).neighbors = [props(ii).neighbors, jj];
            props(ii).dists = [props(ii).dists, dd];
        end
        
        
    end
end


merge = zeros(length(props),1);

for ii = 1:length(props)
    merge(ii, 1) = ii;
    
    if isempty(props(ii).neighbors)
        continue;
    end
    
    ind = ii;
    counter = 1;
    while 1
        [mm, closest_forward] = min(props(ind).dists);
        closest_forward = props(ind).neighbors(closest_forward);
        
        if ~isempty(intersect(merge(ii,:), closest_forward))
            break;
        end
        
        if isNuke(props, nonzeros(unique([merge(ii,:), closest_forward])), 0, MAX_AXIS_VOLUME, MAX_NUKE_DIA, sz )
            counter = counter+1;
            ind = closest_forward;
            merge(ii,counter) = ind;
            
        else
            break;
        end
        
    end
end

merge_ = merge;

counter = 0;
nukes = [];

for ii = size(merge,2):-1:1
    inds = find(sum(~~merge,2)==ii);
    
    while ~isempty(inds)
        ind = find(sum(~~merge,2)==ii, 1, 'first');
        regs = nonzeros(merge(ind,:));
        
        if ~isempty(regs)
            is_circle = 1;
            
            for kk = 1:length(regs)
                regs_regs = nonzeros(merge(regs(kk),:));
                
                if ~isempty(setxor(regs, regs_regs))
                    is_circle = 0;
                    break;
                end
                
            end
            
            if is_circle
                counter = counter+1;
                nukes(counter,1:ii) =  regs;
            
                for kk = 1:length(regs)
                    merge(merge==regs(kk)) = 0;
                end

            else
                counter = counter+1;
                nukes(counter,1) =  regs(1);
                merge(regs(1),:) = 0;
            end
        end
        
        inds = find(sum(~~merge,2)==ii);

    end
end

% setxor(unique(mkL(:)), unique(nukes))

counter = 0;
mkF = 0*mkL;

for ii = 1:size(nukes,1)
    
    inds = nonzeros(unique(nukes(ii,:)));
    [rr,cc,zz] = getBoundingBox3(props, inds, 0, sz);
    
    include_flag = isNuke(props, inds, 0, MAX_AXIS_VOLUME, MAX_NUKE_DIA, sz);
    
    if ~include_flag
        '';
    end
    
    if include_flag
        
        counter = counter + 1;
        
        mkkF = mkF(rr,cc,zz);
        mkkB = mkL(rr,cc,zz);
        
        for iii = 1:length(inds)
            mkkF(mkkB==inds(iii)) = counter;
        end
        
        mkF(rr,cc,zz) = mkkF;
        
    end
end

end




















