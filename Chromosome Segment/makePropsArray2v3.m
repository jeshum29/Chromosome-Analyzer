function [props, status] = makePropsArray(nd, surfsToInclude)


if isfield(nd, 'optiData')
    status = nd.optiData.status;
elseif isfield(nd, 'surf_status')
    status = nd.surf_status;
else
    status = [];
end

if ~exist('surfsToInclude', 'var')
    surfsToInclude = 1:length(nd.sd);
end

numS = length(surfsToInclude);
for ii = 1:numS
    
    r_labels = nd.sd{surfsToInclude(ii)}.neighbors;
    [rr,cc,zz] = getBoundingBox3(nd.ppR, r_labels, 3, size(nd.im));
    
    mkL = nd.sd{surfsToInclude(ii)}.mkL;
    
    [vecs, vals] = calcPrincipleAxes(nd.surfL(rr,cc,zz)==surfsToInclude(ii));
    
    % z-component of the surface's principle axis
    surf_eig_vec_Z(ii)            = abs(vecs(3,1));
    
    % smallest eigenvalue (large if surface is not planar)
    surf_eig_val_minor(ii)        = min(vals);
    
    % intensity along the surface
    surface_mean_intensity(ii)    = mean(nd.im(nd.surfL==surfsToInclude(ii)));
    
    % intensity within the neighboring regions
    region_mean_intensity(ii)     = mean(double(nd.sd{surfsToInclude(ii)}.im(~~mkL)));
    region_std_intensity(ii)      = std(double(nd.sd{surfsToInclude(ii)}.im(~~mkL)));
    
    % region areas
    regs  = nd.sd{surfsToInclude(ii)}.neighbors;
    areas = [nd.ppR(regs).Area];
    
    reg_area(ii)                  = sum(areas);
    reg_ratio(ii)                 = min(areas)/max(areas);
    min_area(ii)                  = min(areas);
    
    
    % approximation to the derivative of the intensity normal to the surface
    reg_sum   = sum(sum(nd.ssd.reg_kymo(surfsToInclude(ii)).kymo, 3),2);
    surf_sum  = sum(sum(nd.ssd.surf_kymo(surfsToInclude(ii)).kymo,3),2);
    surf_pos  = round(sum(surf_sum.*(1:length(surf_sum))')/sum(surf_sum));
    reg_sum   = round(reg_sum(max(1,surf_pos-4):min(end,surf_pos+4)));
    
    reg_dd = (reg_sum(2:end) - reg_sum(1:end-1))./reg_sum(2:end);
    reg_dd(reg_sum(2:end)==0) =  0;
    midd = round(length(reg_dd)/2);
    
    kymo_param(ii) = sum(reg_dd(midd:end)) - sum(reg_dd(1:midd));
    
    
end

% nucleus-wide intensity 
nuke_mean_intensity  = mean(double(nd.im(~~nd.mkAxis(:))));
nuke_std_intensity   = std( double(nd.im(~~nd.mkAxis(:))));

% FUTURE USE
% parameters calculated from already-segmented and traced axes
axis_lengths = 0;
axis_volumes = 0;

if isfield(nd, 'sdata')
    
    trace_cat = [];
    for jj = 1:length(nd.sdata)
        trace_cat = [trace_cat; nd.sdata(jj).trace];
        len(jj) = nd.sdata(jj).length;
    end

    axis_lengths(1:length(len)) = len;

%     [cen, rad] = fitAxesToSphere(trace_cat);
%     nukeRad = (rad(1)*rad(2)*rad(3))^(1/3);

    if max(nd.mkL_final(:))==1
        nd.mkL_final = bwlabeln(nd.mkL_final, 6);
    end

    pp = regionprops(nd.mkL_final, 'Area');
    vol = [pp.Area];
    axis_volumes(1:length(vol)) = vol;

end

surfAreas = [nd.ppS.Area];

props(1,1:numS) = surface_mean_intensity;
props(2,1:numS) = surface_mean_intensity / nuke_mean_intensity;
props(3,1:numS) = surface_mean_intensity ./ region_mean_intensity;
props(4,1:numS) = surf_eig_vec_Z;
props(5,1:numS) = surf_eig_val_minor;
props(6,1:numS) = surfAreas(surfsToInclude);
props(7,1:numS) = reg_area;
props(8,1:numS) = reg_ratio;
props(9,1:numS) = min_area;
props(10,1:numS) = abs(kymo_param);

props(11,1:length(axis_lengths)) = axis_lengths;
props(12,1:length(axis_volumes)) = axis_volumes;














