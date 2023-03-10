function axisData = calcFocusPositionsRAD51(axisData, OFFSET_VECTOR, GENERATE_SPOT_IDS)

% RAD51 focus distribution analysis from OMX images
% Using snakeTrace3v1 axis traces

% Keith Cheveralls
% September 2013


% for 2013_06_24_N2_IF DV data:
%
% manualExcludeList = [21, 23, 26, 27, 77, 83, 125, 127, 138, 139,...
%                      141, 157, 168, 182, 183, 195, 30, 46, 54, 69, 80,...
%                      118, 226, 239, 244, 247, 270, 271, 277, 279];

manualExcludeList = [];

% this is for backwards compatibility
% OFFSET_VECTOR argument used to be only z offset 
% (since x-y alignment was correct until met-2 data on 3/31/16)
if length(OFFSET_VECTOR)==1
    X_CORRECTION_OFFSET = 0;
    Y_CORRECTION_OFFSET = 0;
    Z_CORRECTION_OFFSET = OFFSET_VECTOR;
else
    X_CORRECTION_OFFSET = OFFSET_VECTOR(1);
    Y_CORRECTION_OFFSET = OFFSET_VECTOR(2);
    Z_CORRECTION_OFFSET = OFFSET_VECTOR(3);
end

if ~exist('GENERATE_SPOT_IDS', 'var')
    GENERATE_SPOT_IDS = 0;
    'WARNING: SPOT IDS NOT GENERATED'
end

generateSpotIDsFlag = GENERATE_SPOT_IDS;

spotChannels = [];
if isfield(axisData, 'spots_im1')
    spotChannels = [spotChannels, 1];
end
if isfield(axisData, 'spots_im2')
    spotChannels = [spotChannels, 2];
end
if isfield(axisData, 'spots_im3')
    spotChannels = [spotChannels, 3];
end

if generateSpotIDsFlag
    
    for ii = 1:length(axisData)
        axisData(ii).skipFlag = 0*(1:length(spotChannels));
    end
    
    for kk = 1:length(spotChannels)
        
        spotFieldName     = ['spots_im' num2str(spotChannels(kk))];
        spotFieldNameRedo = [spotFieldName '_redo'];
        
        lastSpotID = 0;
        
        for ii = 1:length(axisData)
            
            % this flag is updated within this loop -- see below
            if axisData(ii).skipFlag(kk)
                continue;
            end
            
            if isfield(axisData(ii), spotFieldNameRedo)
                spots = axisData(ii).(spotFieldNameRedo);
            else
                spots = axisData(ii).(spotFieldName);
            end
            
            numSpots = length(spots);
            
            spotIDs = (lastSpotID + 1):(lastSpotID + numSpots);
            axisData(ii).spotIDs(kk, 1:length(spotIDs)) = spotIDs;
            
            % the entries in this flag vector will be set to zero if the
            % corresponding focus turns out to be closer to a different axis
            axisData(ii).spotAxisAssociationFlag(kk,1:numSpots) = (1:numSpots)*0 + 1;
            
            lastSpotID = lastSpotID + numSpots;
            
            % here we look for axes that came from the same
            % nucleus and assign the same spotIDs to those foci
            % (since they're the same)
            
            for jj = (ii + 1):length(axisData)
                if strcmp(axisData(ii).fname, axisData(jj).fname)
                    
                    spotIDs = axisData(ii).spotIDs(kk,:);
                    axisData(jj).spotIDs(kk, 1:length(spotIDs)) = spotIDs;
                    axisData(jj).skipFlag(kk) = 1;
                end
            end
        end
    end
end

for ii = 1:length(axisData)
    axisData(ii).spotAxisDistances = [];
    axisData(ii).spotAxialPosition = [];
    axisData(ii).totalAxisLength   = [];
    axisData(ii).spotScore         = [];
end

for kk = 1:length(spotChannels)
    
    spotFieldName     = ['spots_im' num2str(spotChannels(kk))];
    spotFieldNameRedo = [spotFieldName '_redo'];
    
    for ii = 1:length(axisData)
        
        if ~isempty(intersect(manualExcludeList, ii))
            continue;
        end
        
        trace = axisData(ii).snakeTrace;
        arcLen = calcArcLength(trace);
        
        if isfield(axisData(ii), spotFieldNameRedo)
            spots = axisData(ii).(spotFieldNameRedo);
        else
            spots = axisData(ii).(spotFieldName);
        end
        
        if isempty(spots)
            axisData(ii).totalAxisLength(kk,1)   = arcLen(end);
        end
        
        for jj = 1:length(spots)
            
            spotPosition = [...
                spots(jj).r(1) - min(axisData(ii).rr) - X_CORRECTION_OFFSET,...
                spots(jj).r(2) - min(axisData(ii).cc) - Y_CORRECTION_OFFSET,...
                spots(jj).r(3) - min(axisData(ii).zz) - Z_CORRECTION_OFFSET,...
                ];
            
            if isempty(trace) || isempty(spotPosition)
                'WARNING in calcFocusDistroRAD51: variable to or r is empty';
            else
                dd = sqrt(sum((trace - repmat(spotPosition, size(trace,1), 1)).^2,2));
                
                [spotAxisDistance,ind] = min(dd);
                
                axisData(ii).totalAxisLength(kk, jj)   = arcLen(end);
                axisData(ii).spotAxisDistances(kk, jj) = spotAxisDistance;
                axisData(ii).spotAxialPosition(kk, jj) = arcLen(ind);
                axisData(ii).spotScore(kk, jj)         = spots(jj).intensity_score;
                
                if length(spotChannels)==1
                    axisData(ii).spotAxisVector(jj, :)   = trace(ind,:) - spotPosition;
                    axisData(ii).spotAbsPosVector(jj, :) = spotPosition;
                else
                    axisData(ii).spotAxisVector(kk, jj, :)   = trace(ind,:) - spotPosition;
                    axisData(ii).spotAbsPosVector(kk, jj, :) = spotPosition;
                end
                
                if 0
                    vizAxisRAD51(axisData, ii);
                    '';
                end
            end
        end
    end
end

return;

% for ii = 1:length(axisData)
%
%     if isempty(axisData(ii).totalAxisLength)
%         continue;
%     end
%
%     for jj = 1:length(axisData(ii).spotAxialPosition)
%
%         spotPosition(ii,jj) = axisData(ii).spotAxialPosition(jj);
%         spotDistance(ii,jj) = axisData(ii).spotAxisDistances(jj);
%         spotScore(ii,jj)    = axisData(ii).spotScore(jj);
%         spotID(ii,jj)       = axisData(ii).spotIDs(jj);
%
%     end
% end

% ------------------------------------------------------------------
%
% ensure that each spot is assigned to only one axis
% (the spotIDs uniquely identify each focus)
%
% This code also appears in plotFcousDistributionv3

% This code is no longer necessary now that spot-axis associations
% are manually supervised using axisViewer - KCC July 2015
% ------------------------------------------------------------------
% 
% if isempty(spots)
%     return;
% end
% 
% uniqueIDs = 1:max(spotID(:));
% 
% for ID = uniqueIDs
%     
%     spotInds = find(spotID(:)==ID);
%     
%     [rows, cols] = ind2sub(size(spotID), spotInds);
%     
%     % if the spot is only associated with one axis, we do nothing
%     if length(spotInds)==1
%         continue;
%     end
%     
%     % otherwise, remove the spot-axis associations for all but
%     % the axis closest to the spot by settings its score to -1
%     
%     dists = spotDistance(spotInds);
%     [min_, min_ind] = min(dists);
%     
%     for jj = 1:length(spotInds)
%         if jj~=min_ind
%             axisData(rows(jj)).spotAxisAssociationFlag(cols(jj)) = 0;
%         end
%     end
% end
% 
