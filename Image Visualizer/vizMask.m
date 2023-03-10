function vizMask(mask, cc, figNum, clfFlag, use_lighting, alpha)

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

if isempty(cc)
    cc = [1,1,1];
end

figure(figNum);

if clfFlag
    clf;
end

 % fv = isosurface(imfilter(double(~~mask),2)>9/10);
 fv = isosurface(double(~~mask));
 p = patch(fv);

 set(p,...
    'facecolor', cc,...
    'edgealpha', 0,...
    'facealpha', alpha);
        
daspect([1 1 1]); 
view(3); 
axis tight; 

axis vis3d;
lighting flat;
material metal;

rotate3d;

if use_lighting
    camlight('headlight')
    camorbit(180,0)
    camlight('headlight')
end


set(gcf, 'color', [0,0,0]);
set(gca, 'Position', [0,0,1,1]);
'';