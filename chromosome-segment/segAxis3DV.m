function nd = segAxis3(nd)

if ~isfield(nd, 'PX_SZ')
    PX_SZ = 64;
    'Warning: no pixel size specified. Using 64nm.'
else
    PX_SZ = nd.PX_SZ;
end

% ============================================================
%
% function to segment chromosomal axes from 3D DV images
%
% Keith Cheveralls
% Dernburg Lab
% June 2013
%
% ============================================================

SEG_CHANNEL    = nd.segParams.SEG_CHANNEL;
SEG_THRESH     = nd.segParams.SEG_THRESH;
SEG_SMOOTHING  = nd.segParams.SEG_SMOOTHING;
SEG_GAMMA      = nd.segParams.SEG_GAMMA;
MIN_SURF_AREA  = nd.segParams.MIN_SURF_AREA;
MIN_REG_AREA   = nd.segParams.MIN_REG_VOL;

% ============================================================--
%
% choose which channel to segment
%

nd.im = eval(['nd.im' num2str(SEG_CHANNEL) ';']);

% ============================================================--

mkNuke = (nd.mk);

sz = size(nd.im);

imDF{1} = imfilter3(nd.im, [1,1,1]*SEG_SMOOTHING/2);
imDF{2} = imfilter3(nd.im, [1,1,1]*SEG_SMOOTHING);

% 
% % ============================================================--
% %
% % make a background mask by Otsu thresholding
% %
% % ============================================================--
% 
% 
% imForMask   = double(nd.im).^SEG_GAMMA;
% SEG_THRESH  = graythresh(autogain(imForMask));
% imForMask   = double(imForMask);
% SEG_THRESH  = SEG_THRESH*max(imForMask(:)) / mean(imForMask(:));
% 
% nd.segParams.SEG_THRESH = SEG_THRESH;
% 
% 
% mkAxisBack = imDF{1} > SEG_THRESH*mean(nd.im(:));

% mkAxisBack = mkAxisBack .* imdilate(mkNuke,ones(5,5,5));

nd = makeAxisBackMask(nd);

% ============================================================--
%
% watershed
%
% ============================================================--

mkAxis = nd.mkAxis;

seeds = imregionalmin(imcomplement(imDF{2}));
ws = watershedManualSeeds(imcomplement(imDF{1}), seeds, 26);

surf = (~ws).*mkAxis;
reg = (~~ws).*mkAxis;

nd.surfRaw = surf;
nd.regRaw = reg;
nd.sk = skel3(double(mkAxis));

% ============================================================-----
%
% preprocess: turn off small regions and turn on small surfaces
%
% ============================================================-----

[reg, regL, ppR, surf, surfL, ppS, sd] = preprocessRegions(reg, surf, MIN_SURF_AREA, MIN_REG_AREA);

nd.surf = surf;
nd.reg = reg;
nd.surfL = surfL;
nd.regL = regL;
nd.sd = sd;
nd.ppR = ppR;
nd.ppS = ppS;


iis = 1:length(ppS);

for ii = iis
    
    r_labels = sd{ii}.neighbors;
    [rr,cc,zz] = getBoundingBox3(ppR, r_labels, 3, sz);
    
    mkk = ...
        (regL(rr,cc,zz)==r_labels(1))*r_labels(1) + ...
        (regL(rr,cc,zz)==r_labels(2))*r_labels(2) + ...
        (surfL(rr,cc,zz)==ii)*0.5;
    
    skk = nd.sk(rr,cc,zz);
    imm = nd.im(rr,cc,zz);
    surf = mkk==1/2;
    
    sd{ii}.mk = ~~mkk;
    sd{ii}.mkL = mkk;
    sd{ii}.im = imm;
    sd{ii}.surf = surf;
    
    mkk = double(~~mkk);
    core = mkk*0;
    delete_exposed_pixels_3(double(~~mkk), core, 26);
    
    shell = double(~~mkk) - core;
    sd{ii}.shell = shell;
    
end


nd.sd = sd;

surfSlicerData = surfSlicer(nd);

nd.ssd = surfSlicerData;

nd.surf_status = iis*0;



