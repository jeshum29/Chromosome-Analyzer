function nd = makeAxisBackMask(nd)

if ~isfield(nd.segParams, 'useGamma')
    nd.segParams.useGamma = 0;
end

if ~isfield(nd.segParams, 'useOtsuThresh')
    nd.segParams.useOtsuThresh = 0;
end

% backwards compatibility
if isfield(nd.segParams, 'GAMMA')
    nd.segParams.SEG_GAMMA = nd.segParams.GAMMA;
end

if ~isfield(nd.segParams, 'SEG_GAMMA')
    nd.segParams.SEG_GAMMA = 1;
end


SEG_CHANNEL    = nd.segParams.SEG_CHANNEL;
SEG_THRESH     = nd.segParams.SEG_THRESH;
SEG_SMOOTHING  = nd.segParams.SEG_SMOOTHING;
SEG_GAMMA      = nd.segParams.SEG_GAMMA;

% ============================================================--
%
% choose which channel to segment
%

im = eval(['nd.im' num2str(SEG_CHANNEL) ';']);
im = double(im);

% ============================================================--

if (nd.segParams.useGamma)
    im = double(im).^SEG_GAMMA;
end

imSmooth = imfilter3(im, [1,1,1]*SEG_SMOOTHING/2);


% ============================================================--
%
% make a background mask by Otsu thresholding
%
% ============================================================--

if nd.segParams.useOtsuThresh
    
    SEG_THRESH = graythresh(autogain(im));
    SEG_THRESH = SEG_THRESH*max(im(:)) / mean(im(:));
    SEG_THRESH = round(100*SEG_THRESH)/100;
    
    nd.segParams.SEG_THRESH = SEG_THRESH;

end

mkAxisBack = imSmooth > SEG_THRESH*mean(im(:));


mkNuke = (nd.mk);
mkAxisBack = mkAxisBack;% .* imdilate(mkNuke,ones(12,12,12));

nd.mkAxis = mkAxisBack;

