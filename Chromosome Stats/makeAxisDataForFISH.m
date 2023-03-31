function axisData = makeAxisMask(nd)

mkL = nd.mkL_final;

axisFlags = nd.mkL_final_includeFlags;

pp = regionprops(mkL, 'BoundingBox');
sz = size(mkL);

goodAxisLabels = find(axisFlags);

axisData = [];

for ii = 1:length(goodAxisLabels)
    
    label = goodAxisLabels(ii);
    
    [rr,cc,zz] = getBoundingBox3(pp, label, 1, sz);
    
    mkk = ~~(mkL(rr,cc,zz)==label) + ~~nd.surfOff(rr,cc,zz);
    
    mkkL = bwlabeln(mkk, 6);
    
    pp_tmp = regionprops(mkkL, 'Area');
    
    aa = max([pp_tmp.Area]);
    
    for jj = 1:length(pp_tmp)
        if pp_tmp(jj).Area~=aa
            mkk(mkkL==jj) = 0;
        end
    end
    
    mkk = imfill(mkk,'holes');
    
    skk = ~~skel3(double(mkk));
    
    skkClean = [];
    try
        skkClean = cleanSkel2(skk);
    end
    
    if isempty(skkClean)
        
        mkk = imerode(mkL(rr,cc,zz)==label, ones(3,3,3));
        skk = ~~skel3(double(mkk));
        try
            skkClean = cleanSkel2(skk);
        end
    end
    
    
    im = nd.im(rr,cc,zz);
    
    thresh = graythresh(autogain(im));
    
    mkkv2       = (im > thresh).*(mkk);
    skkv2       = ~~skel3(double(mkkv2));
    
    skkv2Clean = [];
    try
        skkv2Clean  = cleanSkel2(skkv2);
    end
    
    axisData(ii).axisChannel = nd.segParams.SEG_CHANNEL;
    
    if isfield(nd, 'focusFittingParams');
        axisData(ii).spotChannel = nd.focusFittingParams.channel;
    else
        axisData(ii).spotChannel = [];
    end
    
    axisData(ii).fname  = nd.dir;
    
    axisData(ii).label   = label;
    axisData(ii).sk      = skk;
    axisData(ii).skC     = skkClean;
    axisData(ii).skv2    = skkv2;
    axisData(ii).skv2C   = skkv2Clean;
    axisData(ii).mk      = mkk;
    axisData(ii).mkv2    = mkkv2;
    

    if isfield(nd, 'im1')
        axisData(ii).im1            = nd.im1(rr,cc,zz);
        axisData(ii).im1_globalMean = mean(nd.im1(:));
    end
    
    if isfield(nd, 'im2')
        axisData(ii).im2            = nd.im2(rr,cc,zz);
        axisData(ii).im2_globalMean = mean(nd.im2(:));
    end
    
    if isfield(nd, 'im3')
        axisData(ii).im3            = nd.im3(rr,cc,zz);
        axisData(ii).im3_globalMean = mean(nd.im3(:));
    end
    
    

    if isfield(nd, 'spots_im1')
        axisData(ii).spots_im1      = nd.spots_im1;
    end
    
    if isfield(nd, 'spots_im2')
        axisData(ii).spots_im2      = nd.spots_im2;
    end
    
    if isfield(nd, 'spots_im3')
        axisData(ii).spots_im3      = nd.spots_im3;
    end
    
    
    axisData(ii).rr = rr;
    axisData(ii).cc = cc; 
    axisData(ii).zz = zz;
    
    
end


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

