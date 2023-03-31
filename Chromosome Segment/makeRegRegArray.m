function [regRegArray, surfRegArray, regSurfArray] = makeRegRegArray(regionData, surfaceData, regions, surfaces)

% the ith row of segSegArray lists the segments that neighbor the ith
% segment in the segs array

%regions and surfaces are a list of the regions and surfaces to consider
%(that are not "off")


% create regionData if not supplied

if isempty(regionData)
    
    for ii = 1:length(regions)
        nnbs = [];
        
        for jj = 1:length(surfaces)
            if sum(surfaceData{surfaces(jj)}.neighbors == regions(ii))
                nnbs = [nnbs, surfaces(jj)];
            end
        end
        
        regionData{regions(ii)}.neighbors = nnbs;
        
    end
end


regRegArray = [];

for ii = 1:length(regions)
    
    surfNeighbors = intersect(surfaces, regionData{regions(ii)}.neighbors);
    n = [];
    
    for jj = 1:length(surfNeighbors);
        
        n = [n, surfaceData{surfNeighbors(jj)}.neighbors];
        
    end
    
    n = unique(n);
    n = setdiff(n, regions(ii));
    
    if isempty(n)
        n = 0;
    end
    
    regRegArray(regions(ii), 1:length(n)) = n;
    
end

surfRegArray = [];

for ii = 1:length(surfaces)
    
    counter = 0;
    
    surfRegArray(surfaces(ii),1:length(surfaceData{surfaces(ii)}.neighbors)) = surfaceData{surfaces(ii)}.neighbors;
    
end

regSurfArray = [];

for ii = 1:length(regions)
    
    inds = find(surfRegArray==regions(ii));
    [surfs, yy] = ind2sub(size(surfRegArray), inds);
    
    regSurfArray(regions(ii), 1:length(surfs)) = surfs;
    
    
end
end


