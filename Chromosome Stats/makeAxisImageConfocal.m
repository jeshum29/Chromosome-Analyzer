function imRGB = makeAxisImageFISH(axisData, ind)

fname = axisData(ind).fname;
fname = fname(1: (regexp(fname, 'nukes') - 1));

gonad = load([fname 'gonad.mat']);
gonad = gonad.gonad;



globalMean = [0,0,0];

% 



% if isfield(axisData, 'im1')
%     im1 = axisData(ind).im1;
%     globalMean(1) = axisData(ind).im1_globalMean;
% else
%     'WARNING in vizAxisFISH: no im1 found!'
%     return;
% end
% 
% if isfield(axisData, 'im2')
%     im2 = axisData(ind).im2;
%     globalMean(2) = axisData(ind).im2_globalMean;
% else
%     im2 = 0*im1;
% end
% 
% if isfield(axisData, 'im3')
%     im3 = axisData(ind).im3;
%     globalMean(3) = axisData(ind).im3_globalMean;
% else
%     im3 = 0*im1;
% end

mk = axisData(ind).mk;

im1 = axisData(ind).im;
im2 = 0*im1;
im3 = 0*im2;

gonad.im1_index = [1 2 3];
gonad.im2_index = [];
gonad.im3_index = [];

globalMean = [0,0,0];
globalMean(1) = mean(im1(:));

scale   = [2,2,2];
offset  = [2,2,2];
gamma   = [1,1,1];

% --------------------------------------------------
%
%
% MODIFICATIONS to see the axis in context
%
% 
% offset(1) = offset(1)/4;
% 
% mk = 0*mk + 1;
%
%
% --------------------------------------------------


im1 = scale(1) * (im1 - globalMean(1) * offset(1));
im2 = scale(2) * (im2 - globalMean(2) * offset(2));
im3 = scale(3) * (im3 - globalMean(3) * offset(3));

im1 = double(im1).*mk;
im2 = double(im2).*mk;
im3 = double(im3).*mk;

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


imRGB = cat(4, imr, img, imb);

