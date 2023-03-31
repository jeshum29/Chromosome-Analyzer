function nd = segAxis3(nd, SEG_CHANNEL, SEG_THRESH, SEG_SMOOTHING, MIN_SURF_AREA, MIN_REG_AREA)

if ~exist('SEG_CHANNEL' ,'var')
    SEG_CHANNEL = 1;
end

if ~exist('SEG_SMOOTHING', 'var')
    SEG_SMOOTHING = 1;
end

if ~exist('SEG_THRESH','var')
    SEG_THRESH = 1.5;
end

if ~isfield(nd, 'PX_SZ')
    PX_SZ = 110;
    'Warning: no pixel size specified. Using 110nm.'
else
    PX_SZ = nd.PX_SZ;
end

% ============================================================
%
% function to segment chromosomal axes from 3D confocal images
%
% Keith Cheveralls
% Dernburg Lab
% December 2012
%
% ============================================================

nd.segParams.SEG_CHANNEL = SEG_CHANNEL;
nd.segParams.SEG_THRESH = SEG_THRESH;
nd.segParams.SEG_SMOOTHING = SEG_SMOOTHING;

% ============================================================--
%
% choose which channel to segment
%

nd.im = eval(['nd.im' num2str(SEG_CHANNEL) ';']);

% ============================================================--

mkNuke = (nd.mk);

if PX_SZ > 100

    nd.im = imresize3(nd.im, size(nd.im)*2);
    mkNuke = imresize3(mkNuke, size(mkNuke)*2);

end

sz = size(nd.im);
imDF{1} = imfilter3(nd.im, [1,1,1]*SEG_SMOOTHING/2);
imDF{2} = imfilter3(nd.im, [1,1,1]*SEG_SMOOTHING);

% ============================================================--
%
% make a background mask by mean subtracting
%
% ============================================================--

mkAxisBack = imDF{1} > SEG_THRESH*mean(nd.im(:));
mkAxisBack = mkAxisBack .* imdilate(mkNuke,ones(5,5,5));

% ============================================================--
%
% watershed
%
% ============================================================--

mkAxis = mkAxisBack;

seeds = imregionalmin(imcomplement(imDF{2}));
ws = watershedManualSeeds(imcomplement(imDF{1}), seeds, 26);

surf = (~ws).*mkAxis;
reg = (~~ws).*mkAxis;

nd.mkAxis = mkAxis;
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



