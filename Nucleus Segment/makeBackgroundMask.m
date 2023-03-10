function GONAD = makeBackgroundMask(GONAD)

%----------------------------------------------------------
%
% generate a nucleus mask by threshholding a bandpassed 
% image
%
%----------------------------------------------------------
slash = filesep;

% [dd, ff] = fileparts(dd);
% writeDir = [dd slash ff];
% 
% gonad = load([writeDir slash 'gonad.mat']);
% gonad = gonad.gonad;

dd = GONAD.writeDir;


% use half-sized image for seg to speed things up
imL = double(load3([dd slash 'im' num2str(GONAD.SEG_CHANNEL) '.tif']));

if isfield(GONAD, 'illuminationCompensationProfile')
    
    for ii = 1:size(imL,3)
        imL(:,:,ii) = imL(:,:,ii).*GONAD.illuminationCompensationProfile;
    end
    
end

sz = size(imL);

imL = imresize3(imL, 1/2, 'nearest');

imF_lower = imfilter3(imL, [1,1,1]*GONAD.BANDPASS_HIGH/2);
imF_upper = imfilter3(imL, [1,1,1]*GONAD.BANDPASS_LOW/2);

bandPass = imF_lower - imF_upper;

bandPass = imresize3(bandPass, sz);

save3(autogain16(bandPass), [dd slash], 'mk_bandpass', 0);



% 
% thresher3(bandPass);
% 
% keyboard;

% mkB = bandPass > thresh;
% gonad.bandPass_thresh = thresh;
% 
% save3(uint8(mkB), [dd slash], 'mkB_fullsz', 0);
% save3(autogain(mkB), [dd slash 'mkB'], 'mkB_fullsz_', 1);

% save([dd slash 'bandpass.mat'], 'bandPass');
% save([dd slash 'gonad.mat'], 'gonad');