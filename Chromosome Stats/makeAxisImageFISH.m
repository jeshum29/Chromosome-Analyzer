function imRGB = makeAxisImageFISH(axisData, ind, useAxisMask, useKymoImage)

if ~exist('useAxisMask', 'var')
    useAxisMask = 1;
end
if ~exist('useKymoImage', 'var')
    useKymoImage = 0;
end

imRGB = [];

fname = axisData(ind).fname;
fname = fname(1: (regexp(fname, 'nukes') - 1));

try
    gonad = load([fname 'gonad.mat']);
    gonad = gonad.gonad;
catch
    gonad = axisData(ind).gonad;
end

if ~isfield(gonad, 'im4_scale')
    gonad.im4_scale = 0;
    gonad.im4_off   = 0;
    gonad.im4_gamma = 0;
end

if isstr(gonad.im1_scale)
    gonad.im1_scale  = str2double(gonad.im1_scale);
    gonad.im1_off    = str2double(gonad.im1_off);
    gonad.im1_gamma  = str2double(gonad.im1_gamma);
    gonad.im2_scale  = str2double(gonad.im2_scale);
    gonad.im2_off    = str2double(gonad.im2_off);
    gonad.im2_gamma  = str2double(gonad.im2_gamma);
    gonad.im3_scale  = str2double(gonad.im3_scale);
    gonad.im3_off    = str2double(gonad.im3_off);
    gonad.im3_gamma  = str2double(gonad.im3_gamma);
end

scale   = [gonad.im1_scale, gonad.im2_scale, gonad.im3_scale, gonad.im4_scale];
offset  = [gonad.im1_off,   gonad.im2_off,   gonad.im3_off,   gonad.im4_off];
gamma   = [gonad.im1_gamma, gonad.im2_gamma, gonad.im3_gamma, gonad.im4_gamma];

globalMean = [0,0,0];

if isfield(axisData, 'im1')
    im1 = axisData(ind).im1;
    if useKymoImage
        im1 = real(axisData(ind).kymoImage1);
        im1(im1<0) = 0;
    end
    
    if isfield(axisData, 'im1_globalMean')
        globalMean(1) = axisData(ind).im1_globalMean;
    else
        globalMean(1) = mean(axisData(ind).im1(:));
    end
else
    'WARNING in vizAxisFISH: no im1 found!'
    return;
end

if isfield(axisData, 'im2')
    im2 = axisData(ind).im2;
    if useKymoImage
        im2 = real(axisData(ind).kymoImage2);
        im2(im2<0) = 0;
    end
    if isfield(axisData, 'im2_globalMean')
        globalMean(2) = axisData(ind).im2_globalMean;
    else
        globalMean(2) = mean(axisData(ind).im2(:));
    end
else
    im2 = 0*im1;
end

if isfield(axisData, 'im3')
    im3 = axisData(ind).im3;
    if useKymoImage
        im3 = axisData(ind).kymoImage3;
        im3(im3<0) = 0;
    end
    if isfield(axisData, 'im3_globalMean')
        globalMean(3) = axisData(ind).im3_globalMean;
    else
        globalMean(3) = mean(axisData(ind).im3(:));
    end
else
    im3 = 0*im1;
end

mk = axisData(ind).mk;

if isempty(im1)
    return;
end

% --------------------------------------------------
%
% MODIFICATIONS to see the axis in context
%
% HTP3 in im1
% RAD51 in im2
%
% scale(1:2)   = [1.5, 4];
% offset(1:2)  = [3,   2];
% gamma(1:2)   = [.8, .6];
%
% gonad.im1_index = [2];
% gonad.im2_index = [1];
% gonad.im3_index = [1 2 3];
%
% mk = 0*mk + 1;
% 
% scale = [6,6,0,0];
% offset = [2, 1.7,0,0];
% gamma = [.8, .8,0,0];
%
% --------------------------------------------------
%
% gamma(3) = .8;
% gamma(1) = .6;
% scale(1) = 6;
% offset(1) = 6;

im1 = scale(1) * (im1 - globalMean(1) * offset(1));
im2 = scale(2) * (im2 - globalMean(2) * offset(2));
im3 = scale(3) * (im3 - globalMean(3) * offset(3));

if useAxisMask
    im1 = double(im1).*mk;
    im2 = double(im2).*mk;
    im3 = double(im3).*mk;
end

im1(im1<0) = 0;
im2(im2<0) = 0;
im3(im3<0) = 0;

im1  = autogain(double(im1).^(gamma(1)));
im2  = autogain(double(im2).^(gamma(2)));
im3  = autogain(double(im3).^(gamma(3)));

imr = autogain(0*im1);
img = autogain(0*im1);
imb = autogain(0*im1);


if sum(gonad.im1_index==1)
    imr = imr + im1;
end
if sum(gonad.im1_index==2)
    img = img + im1;
end
if sum(gonad.im1_index==3)
    imb = imb + im1;
end
if sum(gonad.im2_index==1)
    imr = imr + im2;
end
if sum(gonad.im2_index==2)
    img = img + im2;
end
if sum(gonad.im2_index==3)
    imb = imb + im2;
end
if sum(gonad.im3_index==1)
    imr = imr + im3;
end
if sum(gonad.im3_index==2)
    img = img + im3;
end
if sum(gonad.im3_index==3)
    imb = imb + im3;
end


rr =  [0.90196      0.25686      0.25686];
gg =  [ 0.3451      0.98039      0.34902];
bb =  [0.15686      0.58824      0.94118];

rr =  [1 0 0];
gg =  [0 1 0];
bb =  [0 0 1];


imRGB = cat(4,...
    imr*rr(1) + img*gg(1) + imb*bb(1),...
    imr*rr(2) + img*gg(2) + imb*bb(2),...
    imr*rr(3) + img*gg(3) + imb*bb(3));



