
function [include_flag, cen_std] = is_nuke(props, inds, MIN_AXIS_VOLUME, MAX_AXIS_VOLUME, MAX_NUKE_DIA, sz)

a = [];
cen = [];




for ii = 1:length(inds)
    
    a(ii) = props(inds(ii)).Area;

    cen(ii,:) = props(inds(ii)).Centroid;
    
end

   
A_sum = sum(a);

cen_bar = sum( cen .* repmat(a', 1, 3), 1) / A_sum;

cen_bar2 = sum(cen_bar.^2);

cen2_bar = sum(sum( cen .* cen .* repmat(a', 1, 3), 1), 2) / A_sum;


cen_std = sqrt(cen2_bar - cen_bar2);


[rr,cc,zz] = getBoundingBox3(props, inds, 0, sz);

include_flag = ...
    A_sum > MIN_AXIS_VOLUME && ...
    A_sum < MAX_AXIS_VOLUME && ...
    max([rr(end) - rr(1), cc(end) - cc(1), zz(end) - zz(1)]) < MAX_NUKE_DIA && ...
    cen_std < MAX_NUKE_DIA/5;

end