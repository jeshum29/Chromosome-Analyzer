function [surfLabel, sp] = makeSurfsLabel(mk, surf, numReg)

sz = size(mk);
mk = double(mk);

% the connectivity here should be 26! because the watershed is 6-connected

boundary_1 = mk*0;
boundary_2 = mk*0;

label_surface_neighbors(mk, double(surf), boundary_1, boundary_2, 26);
reg = unique(nonzeros([boundary_1(:); boundary_2(:)]));
numRegions = length(reg);

if numRegions < 1000
    offset = 1000;
elseif numRegions >= 1000
    offset = 10000;
end

bb = boundary_1*offset + boundary_2;

inds = unique(nonzeros(bb(:)));
surfLabel = 0*bb;
sp = [];

% relabel surface image

for ii = 1:length(inds)
    surfLabel(bb==inds(ii)) = ii;
    sp{ii}.neighbors = [mod(inds(ii), offset), (inds(ii) - mod(inds(ii),offset))/offset];
end

% in rare cases, two surfaces neighbor the same two regions, 
% which will lead to their receiving the same label

% this code fixes this bug
% 08/10/2012

pp = regionprops(surfLabel,'BoundingBox');

for ii = 1:length(inds)

    [rr,cc,zz] = getBoundingBox3(pp, ii, 5, size(surfLabel));
    surff = surfLabel(rr,cc,zz);
    surffL = bwlabeln(surff==ii, 6);
    num_reg = max(surffL(:));
    
    if num_reg > 1
        
        for jj = 2:num_reg
            
            next_ind = max(surfLabel(:)) + 1;
            
            surff(surffL==jj) = next_ind;
            sp{next_ind}.neighbors = sp{ii}.neighbors;
            
            surfLabel(rr,cc,zz) = surff;
            
        end
    end
end