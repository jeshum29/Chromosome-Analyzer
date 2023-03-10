function snake = snakeTracev1(axisData, ind, opts)

% original snakeTrace algorithm designed to run on axisData
% data structure

% KCC March 2014

disp = 0;

snake = [];

if ~exist('opts', 'var')
    
    opts.pointsPerPixel       = 1;
    
    opts.numIterations        = 30;
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
    opts.pullForceWeight  = .2;
    opts.pullForceOffset  = 2;
    
end

fname = axisData(ind).fname;
fname = fname(1: (regexp(fname, 'nukes') - 1));
gonad = load([fname 'gonad.mat']);
gonad = gonad.gonad;

sk = axisData(ind).skC;
mk = axisData(ind).mk;

% -----------------------------------------------------------------
%
% 
% MODIFICATION MADE 9/27/2013
%
%

mk = imdilate(mk, ones(5,5,5));

%
%
% -----------------------------------------------------------------

try
    im = eval(['axisData(ind).im' num2str(gonad.SEG_CHANNEL)]);
catch
    im = axisData(ind).im;
end

im = double(autogain(im));

% get initial estimate of the curve from the skeleton of the axis volume

skelCoords = [];

try
    skelCoords = getSkelCoordsForSnakeTrace(sk);
catch
    'ERROR in snakeTrace3v1: skeleton is probably circular'
end

if isempty(skelCoords)
    'ERROR in snakeTrace3v1: skeleton is circular or hits the image edge'
    return;
end

% interpolate according to the desired point density
skelCoords = interp1(...
                    1:size(skelCoords,1),...
                    skelCoords,...
                    1:.1:size(skelCoords,1)...
                    );

len = [0; cumsum(sqrt(sum( (skelCoords(2:end,:) - skelCoords(1:end-1,:)).^2, 2)))];

snake = interp1(...
                len,...
                skelCoords,...
                linspace(0, len(end), opts.pointsPerPixel*len(end))...
                );


opts.numPoints  = size(snake,1);
opts.imMean     = mean(im(:));

if opts.numPoints < 3
    'WARNING in snakeTrace3v1: too few points!'
    return;
end

% calculate external energy function from the image of the axis staining
imageEnergy = calcImageEnergy(double(im).*double(~~mk), opts);
imageEnergy = double(autogain(imageEnergy))/255;

% calculate the force by taking the derivative of the energy
imageForce(:,:,:,1) = imageDerivatives3D(imageEnergy, opts.forceSigma, 'x');
imageForce(:,:,:,2) = imageDerivatives3D(imageEnergy, opts.forceSigma, 'y');
imageForce(:,:,:,3) = imageDerivatives3D(imageEnergy, opts.forceSigma, 'z');


% calculate the internal force matrix between the points
internalForceMatrix = calcInternalForceMatrix(opts);


if disp
try    
    imRGB = makeAxisImageFISH(axisData, ind);
catch
    imRGB = makeAxisImageConfocal(axisData, ind);
end

    figure(1);
    clf;
    vizIm(imRGB, 1, []);
    set(gcf, 'Alphamap', (0:64)*.5/64);
    hold on
    
    traceColor = [1/4, 1, 1];
    
    plot3(...
        skelCoords(:,2),...
        skelCoords(:,1),...
        skelCoords(:,3),...
        '-', 'LineWidth', 2,...
        'MarkerSize', 1,...
        'MarkerFaceColor', traceColor,...
        'MarkerEdgeColor', traceColor,...
        'color', traceColor);
end


for ii = 1:opts.numIterations
    
    
    if disp && ii==opts.numIterations
        
        traceColor = [1, 1, 1/4];
        
        plot3(...
            snake(:,2),...
            snake(:,1),...
            snake(:,3),...
            '-', 'LineWidth', 2,...
            'MarkerSize', 1,...
            'MarkerFaceColor', traceColor,...
            'MarkerEdgeColor', traceColor,...
            'color', traceColor);
        
        
        '';
        
    end
    
    makeMovie = 0;
    
    if makeMovie
        
        step = 3;

        angles = [30:step:358, 0:step:28];%, 30:2:358, 0:2:30];

        figure(1);
        clf;
        set(gcf,'Renderer','openGL');
        vizIm(imRGB, 1, []);
        set(gcf, 'Alphamap', (0:64)*.4/64);
        set(gcf, 'color', [0,0,0]);
        
        [a,e] = view;

        hold on
        
        traceColor = [1,1/4,1/4];
        
        plot3(...
            skelCoords(:,2),...
            skelCoords(:,1),...
            skelCoords(:,3),...
            '-', 'LineWidth', 2,...
            'MarkerSize', 1,...
            'MarkerFaceColor', traceColor,...
            'MarkerEdgeColor', traceColor,...
            'color', traceColor);
        
        traceColor = [1/4,1,1/4];
        
        plot3(...
            snake(:,2),...
            snake(:,1),...
            snake(:,3),...
            '-', 'LineWidth', 2,...
            'MarkerSize', 1,...
            'MarkerFaceColor', traceColor,...
            'MarkerEdgeColor', traceColor,...
            'color', traceColor);
        
                
        ff = figure(1);
        view([90, 6]);
        
        set(ff, 'position', [1000, 500, 1200, 900]);
        
        tmp = getframe(ff);
        
        x = size(tmp.cdata,1);
        y = size(tmp.cdata,2);
        
        screenshotR(1:x,1:y,ii) = tmp.cdata(:,:,1);
        screenshotG(1:x,1:y,ii) = tmp.cdata(:,:,2);
        screenshotB(1:x,1:y,ii) = tmp.cdata(:,:,3);
         
        '';
        
    end
    
        
    snake = doSnakeIteration(snake, internalForceMatrix, imageForce, im, opts);
    
    % reinterpolate to enforce constant segment arc length
    if ~mod(ii, opts.reinterpolationRate)
        
        snake_  = interp1(1:size(snake,1), snake, 1:.1:size(snake,1));
        len     = [0; cumsum(sqrt(sum( (snake_(2:end,:) - snake_(1:end-1,:)).^2, 2)))];
        snake   = interp1(len, snake_, linspace(0, len(end), opts.numPoints));
        
    end
    
    
end


'';

if makeMovie
    screenshotRGB = cat(4, screenshotR, screenshotG, screenshotB);    
    save3c(screenshotRGB, 'E:\', 'snakes_movie5', 0);
end

end




