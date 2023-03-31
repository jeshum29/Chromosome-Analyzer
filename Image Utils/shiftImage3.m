function [im] = centerImage(im, positionVec)

%---------------------------------------------------------
%
%---------------------------------------------------------


TDIMS_A = [1,2,3];
TDIMS_B = [1,2,3];
TMAP_B = [];
RESAMP = makeresampler('cubic', 'fill');
F = 0;
sz = size(im);

offset =  positionVec;

T_offset = translateMatrix(offset);

im = tformarray(im,...
    maketform('affine',T_offset),...
    RESAMP, TDIMS_A, TDIMS_B, sz, TMAP_B, F);




