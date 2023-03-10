function vizVolume(sk, cc, figNum, clfFlag, use_lighting, alpha, xOffset, yOffset, zOffset)

if ~exist('figNum', 'var')
    figNum = 1;
    clfFlag = 1;
    use_lighting = 1;
end

if ~exist('use_lighting', 'var')
    use_lighting = 1;
end

if ~exist('alpha', 'var')
    alpha = 1;
end

if ~exist('cc', 'var')
    cc = [1,1,1];
end

if ~exist('xOffset', 'var')
    xOffset = 0;
    yOffset = 0;
    zOffset = 0;
end

if isempty(cc)
    cc = [1,1,1];
end

figure(figNum);

if clfFlag
    clf;
end

cdata = cc;

[xs,ys,zs] = define_cube;

[x,y,z] = ind2sub(size(sk), find(~~sk(:)));

x = x + xOffset;
y = y + yOffset;
z = z + zOffset;

for ii = 1:6
    
    x_tmp = [x' + xs(1,ii); x' + xs(2,ii); x' + xs(3,ii); x' + xs(4,ii)];
    y_tmp = [y' + ys(1,ii); y' + ys(2,ii); y' + ys(3,ii); y' + ys(4,ii)];
    z_tmp = [z' + zs(1,ii); z' + zs(2,ii); z' + zs(3,ii); z' + zs(4,ii)];
    
    p = patch(y_tmp, x_tmp, z_tmp,'b');
    

        set(p,...
            'FaceColor',cdata,...
            'EdgeColor',[0,0,0],...
            'FaceAlpha', alpha,...
            'Facelighting', 'flat',...
            'edgelighting','none',...
            'AmbientStrength', 1,...
            'DiffuseStrength',1,...
            'SpecularStrength', 1);

    '';
    
end

daspect([1 1 1]); 
view(3); 
axis tight; 

axis vis3d;
lighting none; 
 
rotate3d;


if use_lighting
    camlight('headlight')
    camorbit(180,0)
    camlight('headlight')
end

end




function [xs,ys,zs] = define_cube


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


