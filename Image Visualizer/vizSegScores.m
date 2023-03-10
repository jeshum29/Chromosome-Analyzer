function vizInitialSeg(nd, figNum, clfFlag, scores)

if ~exist('figNum', 'var')
    figNum = 1;
    clfFlag = 1;
end

mk = nd.mkAxis;


if isfield(nd, 'optiData')
    status = nd.optiData.status;
else
    status = nd.surf_status;
end

[surfOn, surfOff] = make_axis_mask(nd.surfL, status);

nd.surfOn = surfOn;
nd.surfOff = surfOff;

pp = regionprops(nd.surfL, 'Centroid');

cc = colormap(lines(length(pp)));

% all surfs in different colors
% vizSegColor(nd.surfL, figNum, 1, cc);
vizMask(~~nd.surfOn, [1,0,0], figNum,1,1);
vizMask(~~nd.surfOff, [0,1,0],figNum,0,0);
% mask shells
vizMask(~~nd.regL,[1,1,1]/2, figNum,0,0,1);
% vizMask(mk.*(~nd.surfL), [1,1,1]/2, figNum, 0, 0, .8);

figure(figNum);
set(gcf, 'color', [0,0,0]);
set(gca, 'Position', [0,0,1,1]);

pp_nuke = regionprops(~~nd.mkAxis,'Centroid');

try
    nuke_cen = pp_nuke(1).Centroid;
catch
    nuke_cen = [0,0,0];
end


for ii = (1:length(pp))
    
    offset = 15;
    
    surf_cen = pp(ii).Centroid;
    
    unit_vec = surf_cen - nuke_cen;
    unit_vec = unit_vec ./ sqrt(sum(unit_vec.*unit_vec));
    
    surf_cen = surf_cen + offset * unit_vec;
    
    if 1
        
        text_str = [num2str(ii), ' -- ', num2str(round(scores(ii)*100)/1) '%'];
        cc(ii,:) = [1,0,0]*status(ii) + [0,1,0]*(~status(ii));
        
        text(...
            surf_cen(1), ...
            surf_cen(2), ...
            surf_cen(3), ...
            text_str, 'fontsize', 18, 'color', cc(ii,:));
        
        hold on
        
        plot3(...
            [surf_cen(1), pp(ii).Centroid(1)],...
            [surf_cen(2), pp(ii).Centroid(2)],...
            [surf_cen(3), pp(ii).Centroid(3)],...
            'color', cc(ii,:));
        
    end
    
end


'';

rotate3d;
end

function [surfOn, surfOff] = make_axis_mask(surfL, surf_status)

surfOn = 0*surfL;
surfOff = 0*surfL;

for ii = 1:length(surf_status)
    
    if surf_status(ii)==1
        surfOn(surfL==ii) = 1;
    else
        surfOff(surfL==ii) = 1;
    end
    
end

end

