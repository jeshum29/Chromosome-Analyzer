function nd = displayFinalSeg(nd, figNum)


sz = size(nd.surfL);
surfL = nd.surfL;
reg = nd.reg;
numReg = length(nd.ppR);

cc = nd.colormap;

if ~isfield(nd, 'mkL_final')
    'displayFinalSeg error: no mkL_final field!'
    return;
end

pp_final = regionprops(nd.mkL_final, 'Centroid');

vizSegColor(nd.mkL_final, nd.surfOff, figNum, 1, cc);

pp = regionprops(~~(nd.mkL_final + nd.surfL), 'Centroid');

nuke_cen = pp(1).Centroid;

offset = 40;


FOCI_PRESENT = 0;
if isfield(nd, 'fd')
    if ~isempty(nd.fd)
        FOCI_PRESENT = 1;
    end
end


if FOCI_PRESENT
    
    fd = nd.fd;
    
    for nn = 1:length(fd)
        
        if fd{nn}.sharedFocusIndex
            fd{nn}.maskLabel = fd{fd{nn}.sharedFocusIndex}.maskLabel;
            inds = [fd{nn}.sharedFocusIndex, nn];
        else
            inds = nn;
        end
        
        figure(figNum);
        hold on;
        
        maskLabel = fd{inds(1)}.maskLabel;
        
        chrID = [];
        for jj = 1:length(inds)
            chrID = [chrID, fd{inds(jj)}.chrID];
        end
        
        pos = pp_final(maskLabel).Centroid;
        
        focusPos = fd{nn}.position([2,1,3]);
        
        unit_vec = focusPos - nuke_cen;
        unit_vec = unit_vec ./ sqrt(sum(unit_vec.*unit_vec));
        
        pos = focusPos + offset * unit_vec;
        
        text_str = [];
        for jj = 1:length(chrID)
            ID = chrID(jj);
            switch ID
                case 1
                    text_str = [text_str, ' I'];
                case 2
                    text_str = [text_str, ' II'];
                case 3
                    text_str = [text_str, ' III'];
                case 4
                    text_str = [text_str, ' IV'];
                case 5
                    text_str = [text_str, ' V'];
                case 6
                    text_str = [text_str, ' X'];
                case 0
                    text_str = [text_str, ' del'];
            end
        end
        
        
        text(...
            pos(1),...
            pos(2),...
            pos(3),...
            text_str, 'fontsize', 24, 'color', cc(maskLabel,:));
        
        
        plot3(...
            [pos(1), focusPos(1)],...
            [pos(2), focusPos(2)],...
            [pos(3), focusPos(3)],...
            'color', cc(maskLabel,:));

    end
    
    
    
    
    
    
    
end







