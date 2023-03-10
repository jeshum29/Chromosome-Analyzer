function [fitProps, status] = calcSurfProps(data, dataGlobal)

% -----------------------------------------------------------------------------
%
% Function to calculate population-averaged surface properties from the
% training data set.
% 
% Used by fitSurfScores3 and also optiAxis3SmartFast2v3
%
% Keith Cheveralls
% Dernburg Lab
% March 2013
%
% -----------------------------------------------------------------------------
    


fitProps = [];
status = [];
fitPropsGlobal = [];

% -----------------------------------------------------------------------------
%
% calculate surface properties from the real data set
%
% -----------------------------------------------------------------------------


for ii = 1:length(data)
    
    inds = ~~data(ii).props(1,:);
    
    % this line corrects a bug that occurs when some props(1,:) entries are zero 
    % (which was assumed to be impossible since this is surface mean intensity)
    % KCC 2014-12-10
    inds = 1:size(data(ii).props, 2);
    
    fitProps = [fitProps, [...
        data(ii).props(1,inds);... % mean intensity
        data(ii).props(2,inds);... % mean intensity nucleus-relative
        data(ii).props(3,inds);... % mean intensity region-relative
        ...
        data(ii).props(4,inds);... % surface eig vec z
        data(ii).props(5,inds);... % surface eig val
        data(ii).props(6,inds);... % surface area
        ...
        data(ii).props(7,inds);... % region area  % decent correlation
        data(ii).props(8,inds);... % region ratio  % poor correlation
        data(ii).props(9,inds);... % min area
        ...
        data(ii).props(10,inds);... % kymo param
]];

    if ~isempty(data(ii).status)
        status = [status, data(ii).status(inds)];
    end

end

fitProps(isnan(fitProps)) = 0;

status(status==0) = -1;



% -----------------------------------------------------------------------------
%
% calculate surface properties from the full training data set
%
% -----------------------------------------------------------------------------

for ii = 1:length(dataGlobal)
      
    inds = ~~dataGlobal(ii).props(1,:);
    
    % this line corrects a bug that occurs when some props(1,:) entries are zero 
    % (which was assumed to be impossible since this is surface mean intensity)  
    % KCC 2014-12-10
    inds = 1:size(dataGlobal(ii).props, 2);
    
    fitPropsGlobal = [fitProps, [...
        dataGlobal(ii).props(1,inds);... % mean intensity
        dataGlobal(ii).props(2,inds);... % mean intensity nucleus-relative
        dataGlobal(ii).props(3,inds);... % mean intensity region-relative
        ...
        dataGlobal(ii).props(4,inds);... % surface eig vec z
        dataGlobal(ii).props(5,inds);... % surface eig val
        dataGlobal(ii).props(6,inds);... % surface area
        ...
        dataGlobal(ii).props(7,inds);... % region area
        dataGlobal(ii).props(8,inds);... % region ratio
        dataGlobal(ii).props(9,inds);... % min area
        ...
        dataGlobal(ii).props(10,inds);... % kymo param
]];
    
end

fitPropsGlobal(isnan(fitPropsGlobal)) = 0;

% -----------------------------------------------------------------------------
%
% Here we normalize the surface properties using the training data set mean
%
% -----------------------------------------------------------------------------

fitProps = fitProps ./ repmat(mean(fitPropsGlobal, 2), 1, size(fitProps,2));




