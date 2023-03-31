function vec = rotateVectorBack(vec, tangentVec)

% vec is in the tangent vector's local coordinate system
% which obtained by rotating tangentVec into the 
% x-y-z frame as in getRotatedCenteredImage.m

% this function rotates vec back to the x-y-z frame


[R_1, R_2] = makeRotMatrixFromTangent(tangentVec);

tform_VEC = maketform('affine', R_1*R_2);

vec = tforminv(tform_VEC, vec);
