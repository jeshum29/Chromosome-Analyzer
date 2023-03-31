function [reg, regL, ppR, surf, surfL, ppS, sd] = preprocessRegions(reg, surf, MIN_SURF_AREA, MIN_REG_AREA)

sz = size(reg);

if ~exist('MIN_SURF_AREA', 'var')
    MIN_SURF_AREA = 20;
end
if ~exist('MIN_REG_AREA', 'var')
    MIN_REG_AREA = 200;
end

MIN_ISOLATED_REG_AREA = 1000;


regL = bwlabeln(reg, 26);
ppR = regionprops(regL, 'BoundingBox','Area');

%----------------------------------------------------------
%
% HARD CODED REGION DELETION
%
% first, turn off all really small regions
%
%----------------------------------------------------------

for ii = 1:length(ppR)
    
    if ppR(ii).Area < MIN_REG_AREA
        reg(regL==ii) = 0;
    end
    
end

%----------------------------------------------------------
%
% HARD CODED SURFACE DELETION
%
% note: surface labeling is a bit subtle: always turn on
% surfaces that only border one region, and always turn
% on surfaces bordering more than two regions.
%
%----------------------------------------------------------

surfInts = 0*regL;
count_surface_surface_neighbors(double(regL), double(~~surf), surfInts, 26);

surf(surfInts==1) = 0;
surf(surfInts > 2) = 0;

%----------------------------------------------------------
%
% surf must now be labeled with 6-connectivity to ensure
% that turning off the surface intersections generates
% separate surfaces
%
%----------------------------------------------------------

% bwlabeln(surf, 6);

%----------------------------------------------------------
%
% here we label each surface with an index corresponding to
% the entry in sd that gives the labels of the regions
% each surface borders
%
%----------------------------------------------------------

[surfL, sd] = makeSurfsLabel(regL, surf, length(ppR));

% do bwlabeln(surf, 6) and surfL here have the same labeling - yes
% but a few pixels in bwlabeln(surf, 6) do not appear in surfL. safe to ignore

% from now on, use surfL only, since surfaces that do not border 
% two regions do not appear in surfL, which is desirable

%----------------------------------------------------------
%
% turn on all surfaces that are too small to chop an axis
%
%----------------------------------------------------------

ppS = regionprops(surfL, 'Area','BoundingBox','Centroid');

for ii = 1:length(ppS)
    
    if ppS(ii).Area < MIN_SURF_AREA
        surfL(surfL==ii) = 0;
    end
    
end

%----------------------------------------------------------
%
% turn on all surfaces that separate the same regions at
% two different places. 

% This code does something strange. I think it's written wrong
% KCC 03/01/12
%
%----------------------------------------------------------

% for ii = 1:length(ppS)
%     
%     [rr,cc,zz] = getBoundingBox3(ppS, ii, 3, sz);
%     surff = bwlabeln(surfL(rr,cc,zz)==ii, 6);
%     
%     if max(surff(:)) > 1
%         surfL(surfL==ii) = 0;
%     end
%     
% end

%----------------------------------------------------------
%
% turn of all regions that do not border a surface
%
%----------------------------------------------------------

regList = [];

for ii = 1:length(sd)
    regList = [regList, sd{ii}.neighbors];
end

regList = unique(regList);

for ii = 1:length(ppR)
    
    if isempty(intersect(regList, ii)) && ppR(ii).Area < MIN_ISOLATED_REG_AREA
        reg(regL==ii) = 0;
    end
    
end

%-------------------------------------------------------------------
%
% now relabel surfaces and list the surfaces that border each region
%
%-------------------------------------------------------------------

regL = bwlabeln(reg, 26);
numRegions = max(regL(:));
ppR = regionprops(regL, 'BoundingBox','Area');

[surfL, sd] = makeSurfsLabel(double(regL), double(~~surfL), numRegions);
ppS = regionprops(surfL, 'Area','BoundingBox','Centroid');

surf = ~~surfL;

'';

% 



