function [immR, imm, tangentVecR] = getRotatedCenteredImage(im, mk, crop, positionVec, tangentVec, th_z)

%---------------------------------------------------------
%
% this function rotates a cropped chunk of the image im 
% centered on the pixel at round(positionVec) 
% about the origin given by positionVec (which is not, in general,
% an integer) so that tangentVec lies parallel to the x-axis. 
% 
% mk is an optional mask that 
%---------------------------------------------------------


TDIMS_A = [1,2,3];
TDIMS_B = [1,2,3];
TMAP_B  = [];
RESAMP  = makeresampler('linear', 'fill');
F       = 0;
sz      = size(im);


xx = positionVec(1); 
yy = positionVec(2); 
zz = positionVec(3);

xi = round(xx); 
yi = round(yy); 
zi = round(zz);

rr = (xi - crop(1)) : (xi + crop(1));
cc = (yi - crop(2)) : (yi + crop(2));
zz = (zi - crop(3)) : (zi + crop(3));

cropFlag = ...
    (min(rr) < 1) || ...
    (min(cc) < 1) || ...
    (min(zz) < 1) || ...
    (max(rr) > sz(1)) || ...
    (max(cc) > sz(2)) || ...
    (max(zz) > sz(3));

if cropFlag
    
    'ERROR in getRotatedCenteredImage.m: cropped image is out of bounds'
    return;
    
end

% rr = max(1, xi - crop(1)):min(sz(1), xi + crop(1));
% cc = max(1, yi - crop(2)):min(sz(2), yi + crop(2));
% zz = max(1, zi - crop(3)):min(sz(3), zi + crop(3));
% 

imm = im(rr,cc,zz);

% if ~isempty(mk)
%     mkk = mk(rr,cc,zz);
%     mkk = bwlabeln(mkk,6);
%     label = mkk(crop(1)+1,crop(2)+1,crop(3)+1);
%     
%     if label==0
%         label = mkk(crop(1):crop(1)+2,crop(2):crop(2)+2,crop(3):crop(3)+2);
%         label = median(label(~~label));
%     end
%     
%     if label==0
%         'getRotatedCenteredImage error: no mask density at skeleton position'
%     end
%     
%     imm = imm.*(mkk==label);
% end



if 0
    
    tmp = 0*mk;
    tmptmp = tmp(rr,cc,zz);
    tmptmp(mkk==label) = 1;
    tmp(rr,cc,zz) = tmptmp;

    vizVolume2(~~mk,0,[1,1,1],1,1,1/10);
    
    loc = 0*mk;
    loc(xi,yi,zi) = 1;
    vizVolume2(tmp,0,[0,1,0],1,0,1/4);
    vizVolume2(loc, 0, [1,0,0],1,0,1);
    
    '';
end


center = [ mean(rr), mean(cc), mean(zz) ];


offset = center - positionVec;

T_offset = translateMatrix(offset);

imm = tformarray(...
        imm,...
        maketform('affine',T_offset),...
        RESAMP,...
        TDIMS_A,...
        TDIMS_B,...
        size(imm),...
        TMAP_B,...
        F);

image_center = (size(imm)+1)/2;

T_F = translateMatrix(-image_center);
T_R = translateMatrix( image_center);

if exist('th_z', 'var')
    [R_1, R_2]   = makeRotMatrixFromTangent(tangentVec, th_z);
else
    [R_1, R_2]   = makeRotMatrixFromTangent(tangentVec);
end

tform_IMAGE  = maketform('affine', T_F*R_1*R_2*T_R);
tform_VEC    = maketform('affine', R_1*R_2);

immR = tformarray(...
        imm,...
        tform_IMAGE,...
        RESAMP,...
        TDIMS_A,...
        TDIMS_B,...
        size(imm),...
        TMAP_B,...
        F);

tangentVecR = tformfwd(tform_VEC, tangentVec);


% [THETA_Z, THETA_Y];

'';


