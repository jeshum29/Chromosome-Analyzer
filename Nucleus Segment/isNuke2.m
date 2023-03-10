function is_nuke = isNuke2(area, dia, mkk, MIN_NUKE_VOLUME, MAX_NUKE_VOLUME, MAX_NUKE_DIA, MAX_BORDER_AREA)

    
if ~exist('MAX_BORDER_AREA', 'var')
    MAX_BORDER_AREA = 50;
end

    is_nuke = ...
        area > MIN_NUKE_VOLUME && ...
        area < MAX_NUKE_VOLUME && ...
        dia  < MAX_NUKE_DIA && ...
        (sum(sum(sum(  mkk( :,       [1,end], :       )))) < MAX_BORDER_AREA) &&...
        (sum(sum(sum(  mkk( [1,end], :,       :       )))) < MAX_BORDER_AREA) &&...
        (sum(sum(sum(  mkk( :,       :,       [1,end] )))) < MAX_BORDER_AREA);
    
    
    
    