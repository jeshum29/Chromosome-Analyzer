function [R_1, R_2] = makeRotMatrixFromTangent(tangentVec, th_z)

tangentVec = tangentVec / sqrt(sum(tangentVec.^2));

THETA_Z = atan((tangentVec(2)/tangentVec(1)));

if isnan(THETA_Z)
    THETA_Z = 0;
end

if tangentVec(1) < 0
   % THETA_Z = pi + THETA_Z;
end

if exist('th_z','var')
    THETA_Z = th_z;
end

R_1 = rotMatrix_Z(THETA_Z);

vec_tmp = tformfwd(maketform('affine', R_1), tangentVec);

THETA_Y = atan((vec_tmp(3) / vec_tmp(1)));

if isnan(THETA_Y)
    THETA_Y = 0;
end

if vec_tmp(1) < 0
   % THETA_Y = pi + THETA_Y;
end

R_2 = rotMatrix_Y(THETA_Y);




'';