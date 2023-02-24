
function [surfOn, surfOff] = make_axis_mask(surfL, surf_status)

surfOn = 0*surfL;
surfOff = 0*surfL;

for ii = 1:length(surf_status)
    
    if surf_status(ii)==10
        surfOn = surfOn + (surfL==ii)*2;
    elseif surf_status(ii)==1
        surfOn = surfOn + (surfL==ii);
    else
        surfOff = surfOff + (surfL==ii);
    end
end

end