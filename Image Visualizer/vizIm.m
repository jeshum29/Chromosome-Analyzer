function vizMask(im, figNum, imAlpha, xdata, ydata, zdata)

sz = size(im);


if ~exist('figNum', 'var')
    figNum = 1;
end

if ~exist('imAlpha', 'var')
    imAlpha = [];
end


figure(figNum);

if ~isempty(imAlpha)
    
    if exist('xdata', 'var')
        
        vol3d('Cdata', im, 'texture', '3D', 'Alpha', imAlpha,...
            'XData', xdata, 'YData', ydata, 'ZData', zdata);
    
    else
        
        vol3d('CData', im,'texture','3D','Alpha', imAlpha);
    
    end
else
    vol3d('CData', im,'texture','3D');
end


daspect([1 1 1]);
view(3);
axis tight;
%  camlight left;
%  camlight right;
axis vis3d;
%  lighting none; %without this, no intensity mapping
lighting flat;
material dull;

rotate3d;
set(gcf, 'color', [0,0,0]);
set(gca, 'Position', [0,0,1,1]);

'';