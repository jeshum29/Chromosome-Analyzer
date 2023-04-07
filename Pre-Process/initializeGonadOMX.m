function gonad = initializeGonadOMX(dd, NUM_CHAN, SEG_CHANNEL, OVERWRITE_FLAG)

%-------------------------------------------------------------------------------------
%
%  settings for segmenting confocal images of axis markers
%
%
%  Keith Cheveralls
%  Dernburg Lab
%  July 2013
%
%-------------------------------------------------------------------------------------

if ~exist('OVERWRITE_FLAG', 'var')
    OVERWRITE_FLAG = 0;
end

% here we load images with resize in z to isotropic resolution

gonad.NUM_CHAN = NUM_CHAN;

gonad.isConfocal = 0;

slash = filesep;

% [dd, ff] = fileparts(dd);
% writeDir = [dd slash ff];
% 
% mkdir(writeDir);

writeDir = dd;

if OVERWRITE_FLAG
    
else
    
if exist([writeDir slash 'gonad.mat'], 'file')
    'WARNING: gonad.mat already exists! Aborting mission.'
      return;
end

end


%-------------------------------------------------------------------------------------
%
% CHANNEL TO USE FOR SEGMENTATION
%

gonad.SEG_CHANNEL = SEG_CHANNEL;

%
%-------------------------------------------------------------------------------------

gonad.dd = dd;
gonad.writeDir = writeDir;

% RESAMPLED PIXEL SIZE
gonad.PHYSICAL_PX_SZ = [1,1,1]*79.2/2; %in nm


% nucleus seg settings

gonad.BANDPASS_LOW   = 16;
gonad.BANDPASS_HIGH  = 4.8;

gonad.HIGH_THRESH    = 23000;
gonad.LOW_THRESH     = 23000;

gonad.MIN_AXIS_VOLUME = 1000;
gonad.MAX_AXIS_VOLUME = 120000;

gonad.MIN_NUKE_VOLUME = 35000;
gonad.MAX_NUKE_VOLUME = 120000;

gonad.MAX_BORDER_AREA = 100;

gonad.MAX_NUKE_DIA    = 170;
gonad.MIN_NUKE_DIA    = 80;

gonad.ASPECT_RATIO    = 1;

% segmentation settings

gonad.SEG_SMOOTHING  = 3;
gonad.SEG_THRESH     = 1.7;
gonad.SEG_GAMMA      = .9;

gonad.MIN_REG_VOL    = 50;
gonad.MIN_SURF_AREA  = 5;

gonad.useGamma         = 1;
gonad.useOtsuThresh    = 1;
        
% Display settings

gonad.im1_scale   = '6';
gonad.im2_scale   = '6';
gonad.im3_scale   = '3';
gonad.im1_off     = '2';
gonad.im2_off     = '2';
gonad.im3_off     = '3';
gonad.im1_gamma   = '.8';
gonad.im2_gamma   = '.6';
gonad.im3_gamma   = '.8';
        

'CRITICAL ALERT: Overwriting gonad.mat'
save([writeDir slash 'gonad.mat'], 'gonad');



