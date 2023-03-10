function vizVolume(varargin)

[sk, im, colorMap, edgeColor, figNum, clfFlag, alpha] = parseInput(varargin);

useIntensityImage = 1;

if unique(im(:))==0
    useIntensityImage = 0;
end

if figNum
    ff = figure(figNum);
    if clfFlag
        clf;
    end
end

% set(ff, 'color', [0,0,0]);

if ~useIntensityImage
    cdata = colorMap;
else
    cdata = define_cdata(colorMap, im, sk);
end


[xs,ys,zs] = define_cube;

[x,y,z] = ind2sub(size(sk), find(~~sk(:)));


for ii = 1:6
    
    x_tmp = [x' + xs(1,ii); x' + xs(2,ii); x' + xs(3,ii); x' + xs(4,ii)];
    y_tmp = [y' + ys(1,ii); y' + ys(2,ii); y' + ys(3,ii); y' + ys(4,ii)];
    z_tmp = [z' + zs(1,ii); z' + zs(2,ii); z' + zs(3,ii); z' + zs(4,ii)];
    
    p = patch(y_tmp, x_tmp, z_tmp,'b');
    
    if useIntensityImage
        set(p,...
            'FaceColor','flat',...
            'CData',cdata,...
            'FaceAlpha', alpha,...
            'EdgeColor', edgeColor,...
            'Facelighting', 'flat', ...
            'edgelighting','none',...
            'AmbientStrength', 0,...
            'DiffuseStrength',0);
    else
        set(p,...
            'FaceColor',cdata,...
            'EdgeColor',edgeColor,...
            'FaceAlpha', alpha,...
            'Facelighting', 'flat',...
            'edgelighting','none',...
            'AmbientStrength', 1,...
            'DiffuseStrength',1,...
            'SpecularStrength', 1);
    end
    '';
    
end

daspect([1 1 1]); 

% view(3); 
% axis tight; 
% camlight left; 
% camlight right; 

axis vis3d;
 lighting none; %without this, no intensity mapping
 %lighting flat;
%  material dull;

rotate3d;

end


function cdata = define_cdata(colorMap, im, sk)


if length(size(im))==3
    
    im = autogain(im);
    intensities = double(im(~~sk));
    clear cdata;
    
    if size(colorMap,1) < 256
        
        colorMap = 0:1/255:1;
        colorMap = [colorMap',colorMap',colorMap'];
        
    end
    
    
    if ~isempty(colorMap)
        
        cdata(:,:,1) = colorMap(intensities+1,1);
        cdata(:,:,2) = colorMap(intensities+1,2);
        cdata(:,:,3) = colorMap(intensities+1,3);
        
    else
        
        cdata(:,:,1) = intensities'/255;
        cdata(:,:,2) = intensities'/255;
        cdata(:,:,3) = intensities'/255;
        
    end
    
else
    
    imR = (im(:,:,:,1)); imG = (im(:,:,:,2)); imB = (im(:,:,:,3));
    
    cdata(:,:,1) = (double( imR(~~sk))'/255);
    cdata(:,:,2) = (double(imG(~~sk))'/255);
    cdata(:,:,3) = (double(imB(~~sk))'/255);
    
    cdata = double(autogain(cdata))/255;
    cdata(cdata==0) = 1;
    
    
end


end

function [sk, im, cc, edgeColor, figNum, clfFlag, alpha] = parseInput(varargin)

varargin = varargin{1};
sk = varargin{1};

L = length(varargin);

im = 0;
cc = [1,1,1];
edgeColor = [0,0,0];
figNum = 1;
clfFlag = 1;
alpha = 1/4;

switch L
    
    case 1
        return;
        
    case 2
        
        if length(varargin{2})==3
            cc = varargin{2};
        else
            im = varargin{2};
            
        end
        
    case 3
        
        im = varargin{2};
        cc = varargin{3};
        
    case 4
        im = varargin{2};
        cc = varargin{3};
        edgeColor = varargin{4};
        
    case 5
        im = varargin{2};
        cc = varargin{3};
        clfFlag = varargin{5};
        figNum = varargin{4};
        
    case 6
        im = varargin{2};
        cc = varargin{3};
        alpha = varargin{6};
        figNum = varargin{4};
        clfFlag = varargin{5};
    case 7
        
        im = varargin{2};
        cc = varargin{3};
        edgeColor = varargin{7};
        figNum = varargin{4};
        clfFlag = varargin{5};
        alpha = varargin{6};
        
end

end


function [xs,ys,zs] = define_cube;


xs = [...
    0 1 1 0 0 0;...
    1 1 0 0 1 1;...
    1 1 0 0 1 1;...
    0 1 1 0 0 0]...
    - 0.5;

ys = [...
    0 0 1 1 0 0;...
    0 1 1 0 0 0;...
    0 1 1 0 1 1;...
    0 0 1 1 1 1]...
    - 0.5;

zs = [...
    0 0 0 0 0 1;...
    0 0 0 0 0 1;...
    1 1 1 1 0 1;...
    1 1 1 1 0 1]...
    - 0.5;

end


