function nd = makeAxisMask(nd)


if isfield(nd,'optiData')
    [surfOn, surfOff] = make_axis_mask(nd.surfL, nd.optiData.status);
else
    [surfOn, surfOff] = make_axis_mask(nd.surfL, nd.surf_status);
end

nd.surfOn = surfOn;
nd.surfOff = surfOff;

mk = nd.reg + nd.surfOff;
mk = ~~mk;

% this should remove holes in the axis mask
mk = imfill(mk,'holes');

pp = regionprops(bwlabeln(~mk, 6),'PixelIdxList');

% just in case, remove holes again
for ii = 1:length(pp)
    if length(pp(ii).PixelIdxList) < 1000
        mk(pp(ii).PixelIdxList) = 1;
    end
end

% label the mask
[mkL, num] = bwlabeln(mk, 6);

pp = regionprops(mkL, 'BoundingBox');
sz = size(mk);

sk = 0*mk;
sk_c = 0*mk;

% create a correspondingly labeled skeleton
for ii = 1:num
    
    [rr,cc,zz] = getBoundingBox3(pp, ii, 1, sz);
    
    % make a skeleton
    mkk = mkL(rr,cc,zz)==ii;
    skk = ~~skel3(double(mkk));
    
    % if cleanSkel fails, there's probably a terminal loop or other weird topology, so try dilating the mask
    try 
        skk_c = cleanSkel2(skk);
    catch
        mkk = imdilate(mkL(rr,cc,zz)==ii, ones(3,3,3));
        skk = ~~skel3(double(mkk));
        try
            skk_c = cleanSkel2(skk);
        catch
            skk_c = [];
        end
    end
    
    sk_tmp = sk(rr,cc,zz);
    sk_tmp(~~skk) = ii;
    sk(rr,cc,zz) = sk_tmp;
    
    if ~isempty(skk_c)
        
        sk_tmp = sk_c(rr,cc,zz);
        sk_tmp(~~skk_c) = ii;
        sk_c(rr,cc,zz) = sk_tmp;
    
    end
    
end

nd.mkL_final = mkL;
nd.skL_final = sk;
nd.skL_final_c = sk_c;

try
    nd = rmfield(nd, 'sdata');
end

% now skeleton properties for each axis skeleton
pp = regionprops(nd.skL_final, 'Area','BoundingBox');
[aa, ord] = sort([pp.Area], 'descend');

for ii = 1:min(6,length(ord))
        
    ind = ord(ii);
    
    [rr,cc,zz] = getBoundingBox3(pp,ind,4,sz);
    skk = nd.skL_final(rr,cc,zz)==ind;
    skk_c = nd.skL_final_c(rr,cc,zz)==ind;

    sd.label = ind;
    sd.sk = skk;
    sd.sk_c = skk_c;
    sd.rr = rr; sd.cc = cc; sd.zz = zz;
    sd.length = sum(skk_c(:));

    nd.sdata(ii) = sd;
        
end

nd.sprops = pp;


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

