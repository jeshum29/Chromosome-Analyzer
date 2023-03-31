function opts = defineOpts

% -----------------------------------------------------------------
%
% 
% Original parameters
%
%

opts.pointsPerPixel       = 1;

opts.numIterations        = 80;
opts.reinterpolationRate  = 10;
opts.dt                   = 1;
opts.ds                   = 1;

% extension stiffness
opts.alpha                = .1;

% bending stiffness
opts.beta                 = .03;

% smoothing widths
opts.forceSigma           = 1/2;
opts.energySigma          = 1/2;

% image weights
% positive weight attracts snake to local maxima
opts.imageWeight      =   5;
opts.edgeWeight       =   0;
opts.lineWeight       =   1;

% weights for the pulling forces at snake ends
opts.pullForceWeight  = .1;
opts.pullForceOffset  =  4;


%
%
%
% -----------------------------------------------------------------




% -----------------------------------------------------------------
%
% 
% Exploratory changes made 9/27/2013
%
%

% weights for the pulling forces at snake ends

% .1 looks about right
% .3 looks a bit too extended
% .03 is only slightly extended past the skeleton ends

opts.pullForceWeight  = .1;

% 4 looks about right
% 2 does not look much different

opts.pullForceOffset  =  4;

%
%
%
% -----------------------------------------------------------------

