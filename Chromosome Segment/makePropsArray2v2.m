function props = makePropsArray(nd)

numS = length(nd.sd);

for ii = 1:numS
    
    r_labels = nd.sd{ii}.neighbors;
    [rr,cc,zz] = getBoundingBox3(nd.ppR, r_labels, 3, size(nd.im));
    
    mkL = nd.sd{ii}.mkL;
    
    
    % z-component of the surface's orientation
    [surf_eig_vecs] = calcPrincipleAxes(nd.surfL(rr,cc,zz)==ii);
    surf_eig_vecs_array(:,:,ii) = surf_eig_vecs;
    
    % intensity on along the surface and within the region
    surface_mean_intensities(ii) = mean(nd.im(nd.surfL==ii));
    reg_Ibar(ii) = mean(double(nd.sd{ii}.im(~~mkL)));
    reg_Istd(ii) = std(double(nd.sd{ii}.im(~~mkL)));
    
    % region areas
    regs = nd.sd{ii}.neighbors;
    areas = [nd.ppR(regs).Area];
    reg_areas(ii) = sum(areas);
    reg_ratio(ii) = min(areas)/max(areas);
    min_area(ii) = min(areas);
    
    % approximation to the derivative of the intensity normal to the surface
    reg_sum = sum(sum(nd.ssd.reg_kymo(ii).kymo,3),2);
    surf_sum = sum(sum(nd.ssd.surf_kymo(ii).kymo,3),2);
    [mx, surf_pos] = max(surf_sum);
    reg_sum = reg_sum(max(1,surf_pos-4):min(end,surf_pos+4));
    reg_dd = (reg_sum(2:end) - reg_sum(1:end-1))./reg_sum(2:end);
    reg_dd(reg_sum(2:end)==0) = 0;
    midd = round(length(reg_dd)/2);
    kymo_param(ii) = sum(reg_dd(midd:end)) - sum(reg_dd(1:midd));
    
    
end

% nucleus-wide intensity measures
Ibar = mean(double(nd.im(~~nd.mkAxis(:))));
Istd = std(double(nd.im(~~nd.mkAxis(:))));
intensity_scores = surface_mean_intensities / Ibar;
intensity_scores_reg = surface_mean_intensities / reg_Ibar;

% nucleus volume and sum of surface areas
surfNum = numS;
nukeVol = sum([nd.ppR.Area]);
surfArea = sum([nd.ppS.Area]);

% parameters calculated from already-segmented and traced axes

nukeRad = 0;
axis_lengths = 0;
axis_volumes = 0;

if isfield(nd, 'sdata')
    
    trace_cat = [];
    for jj = 1:length(nd.sdata)
        trace_cat = [trace_cat; nd.sdata(jj).trace];
        len(jj) = nd.sdata(jj).length;
    end

    axis_lengths(1:length(len)) = len;

    [cen, rad] = fitAxesToSphere(trace_cat);
    nukeRad = (rad(1)*rad(2)*rad(3))^(1/3);

    if max(nd.mkL_final(:))==1
        nd.mkL_final = bwlabeln(nd.mkL_final, 6);
    end

    pp = regionprops(nd.mkL_final, 'Area');
    vol = [pp.Area];
    axis_volumes(1:length(vol)) = vol;

end

props(1,1:numS) = nd.surf_status;
props(2,1:numS) = surface_mean_intensities;
props(3,1:numS) = intensity_scores;
props(4,1:numS) = intensity_scores_reg;
props(5,1:numS) = abs(surf_eig_vecs_array(3,1,:));
props(6,1:numS) = [nd.ppS.Area];
props(7,1:numS) = reg_areas;

props(8,:) = surfNum;
props(9,:) = nukeVol;
props(10,:) = surfArea;
props(11,:) = nukeRad;

props(12,1:length(axis_lengths)) = axis_lengths;
props(13,1:length(axis_volumes)) = axis_volumes;

props(14,1:numS) = kymo_param;
props(15,1:numS) = reg_ratio;
props(16,1:numS) = min_area;


