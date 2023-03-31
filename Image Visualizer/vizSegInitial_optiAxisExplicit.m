function vizInitialSeg(nd, figNum, clfFlag, always_off, always_on)

if ~exist('figNum', 'var')
    figNum = 1;
    clfFlag = 1;
end

mk = nd.mkAxis;

use_lighting = 1;



pp = regionprops(nd.surfL, 'Centroid');

cc = colormap(lines(length(pp)));

[surfOn, surfOff, surfPossible] = make_axis_mask(nd.surfL, always_on, always_off);

vizMask(~~surfOn, [1,0,0], figNum,1,1);
vizMask(~~surfOff, [0,1,0],figNum,0,0);
vizMask(~~surfPossible, [1,1,0],figNum,0,0);

vizMask(~~nd.regL,[1,1,1]/2, figNum,0,0,1);


show_labels = 1;

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
                
                if show_labels
                                    
                    text_str = num2str(ii);
                    cc(ii,:) = [1,0,0]*nd.surf_status(ii) + [0,1,0]*(~nd.surf_status(ii));
                    
                text(...
                    surf_cen(1), ...
                    surf_cen(2), ...
                    surf_cen(3), ...
                    text_str, 'fontsize', 18, 'color', cc(ii,:));
                
                hold on
                
%                 plot3(...
%                     [surf_cen(1), pp(ii).Centroid(1)],...
%                     [surf_cen(2), pp(ii).Centroid(2)],...
%                     [surf_cen(3), pp(ii).Centroid(3)],...
%                     'color', cc(ii,:));
                
                end
                
            end
            
            
'';

rotate3d;



function [surfOn, surfOff, surfPossible] = make_axis_mask(surfL, always_on, always_off)

surfOn = 0*surfL;
surfOff = 0*surfL;
surfPossible = 0*surfL;

for ii = 1:length(always_on)
    
    if always_on(ii)==1
        surfOn(surfL==ii) = 1;
    elseif always_off(ii)==1
        surfOff(surfL==ii) = 1;
    else
        surfPossible(surfL==ii) = 1;
    end
    
end

end

end



