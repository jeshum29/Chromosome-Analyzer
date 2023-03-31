function batch_settings(dds)

% 488: RAD-51
% 561: HTP-3

% channel 1 - 488
% channel 2 - cy3


% --------------------------------------------------------
% 
% IMPORTANT NOTE
%
% The illumination compensation profile is different for this data.
% The multiplier, which has been hardcoded to 3 in all other 
% 2015 and 2014 data, is now set to 5 explicitly below.
%
illuminationCompensationProfileIntensity = 5;
%
% KCC June 2015
%
% --------------------------------------------------------

for ii = 1:length(dds)
    
    if ~isdir(dds(ii).name)
        continue;
    end
    
    gonad = load([dds(ii).name filesep 'gonad.mat']);
    gonad = gonad.gonad;
     
   % gonad.LOW_THRESH      = 24000;
   % gonad.MAX_AXIS_VOLUME = 170000;
   % gonad.MAX_NUKE_VOLUME = 170000;
   %  gonad.MAX_NUKE_DIA    = 170;
   %  gonad.MAX_BORDER_AREA = 500;
    
    gonad.im1_scale = 3;
    gonad.im2_scale = 3;
    gonad.im3_scale = 3;
    gonad.im4_scale = 0;
    
    gonad.im1_off   = 3;
    gonad.im2_off   = 3;
    gonad.im3_off   = 5;
    gonad.im4_off   = 0;
    
    gonad.im1_gamma = .7;
    gonad.im2_gamma = .7;
    gonad.im3_gamma = 1;
    gonad.im4_gamma = 0;
    
    gonad.im1_index = [  1   ];
    gonad.im2_index = [  1 2 3   ];
    gonad.im3_index = [  ];
    gonad.im4_index = [];
    
    gonad.chan1_tag = 'RAD51';
    gonad.chan2_tag = 'HTP3';
    gonad.chan3_tag = 'none';
    gonad.chan4_tag = 'none';
    
    % segmentation settings
    gonad.SEG_CHANNEL      = 2;
    gonad.SEG_SMOOTHING    = 3;
    gonad.SEG_GAMMA        = .8;
    
    % SEG_THRESH is not used if useOtsuThresh is enabled
    gonad.SEG_THRESH       = 2;
    gonad.useGamma         = 1;    
    gonad.useOtsuThresh    = 1;
    
    gonad.MIN_REG_VOL      = 25;
    gonad.MIN_SURF_AREA    = 1;
    
    gonad.labelAxesInVolumeOrder = 1;
    
    % focus-finding settings
    % There exist RAD-51 foci only
    fATable = zeros(6,3);
    gonad.focusAxisTable    = fATable;
    gonad.focusChannels     = [0,   1];
    gonad.focusSmoothing    = [0,   1];
    gonad.focusThresh       = [0,  10];
    
    % create an empirical illumination profile to
    % approximately compensate for the inhomogenous
    % illumination by multiplying the intensity image
    
    if ~isfield(gonad, 'illuminationCompensationProfile')
        kernel = fspecial('gaussian', 1024, 800);
        kernel = kernel/max(kernel(:));
        
        kernel = imcomplement(kernel) * illuminationCompensationProfileIntensity + 1;
        
        gonad.illuminationCompensationProfile = kernel;
    end
    
    try
        save( [gonad.writeDir filesep 'gonad.mat'], 'gonad');
    catch
        'WARNING: Failed to save gonad.mat'
    end
    
    
    
end


