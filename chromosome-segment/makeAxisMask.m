function nd = makeAxisMask(nd, labelAxesInVolumeOrder)

%-------------------------------------------------------------------
%
% Function to create a final labeled axis mask
% from a list of surface statuses.
%
% ***IMPORTANT***
% The manual axis-selection part of nukeViewer and downstream
% analysis of selected axes relies on this function always
% generating the same mask from the same status vector.
%
% Keith Cheveralls
% July 2013
%
% Edited December 2014 to allow labeling of axes in volume order
% (to ensure that real axes are selectable in nukeViewer)
%--------------------------------------------------------------------

if ~exist('labelAxesInVolumeOrder', 'var')
    labelAxesInVolumeOrder = 0;
end

regL = nd.regL;
numReg = length(nd.ppR);

if isfield(nd, 'optiData')
    status = nd.optiData.status;
else
    status = nd.surf_status;
end

[surfOn, surfOff] = make_axis_mask(nd.surfL, status);

nd.surfOn = surfOn;
nd.surfOff = surfOff;

mk = nd.reg + nd.surfOff;
mk = ~~mk;

% label the mask
[mkL, num] = bwlabeln(mk, 6);

regRegArray = makeRegRegArray([], nd.sd, 1:numReg, find(~status));
groups = makeSubGroups(regRegArray, 1:numReg);

% relabel the mask according to the group index
mkl_final = 0*mkL;

for ii = 1:length(groups)
    for jj = 1:length(groups(ii).regs)
        mkl_final(regL==groups(ii).regs(jj)) = ii;
    end
end


if labelAxesInVolumeOrder
    
   mkl_final_ = mkl_final*0;
   
   props = regionprops(mkl_final, 'Area');
   
   [ss,ord] = sort([props.Area], 'descend');
   
   for ii = 1:length(ord)
       mkl_final_(mkl_final==ord(ii)) = ii;
   end
   
   mkl_final = mkl_final_;
   
end

nd.mkL_final = mkl_final;
nd.groups_final = groups;

cc = colormap(lines(max(nd.mkL_final(:))));
nd.colormap = cc;

if isfield(nd, 'fd')
    fd = nd.fd;
    groups = nd.groups_final;
    
    for nn = 1:length(fd)
        for ii = 1:length(groups)
            if ~isempty(intersect(groups(ii).regs, fd{nn}.region))
                fd{nn}.regs = groups(ii).regs;
                fd{nn}.maskLabel = ii;
                break;
            end
        end
    end
    
    nd.fd = fd;
end




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
