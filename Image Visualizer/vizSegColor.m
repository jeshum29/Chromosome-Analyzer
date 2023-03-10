function vizSegColor(mk, surfs, figNum, clfFlag, cc)

if ~exist('figNum', 'var')
    figNum = 1;
    clfFlag = 1;
end

% mk = nd.reg + nd.surfOff;
% mk = bwlabeln(mk, 6);

if clfFlag
    try
    close(figNum);
    end
    figure(figNum);
end

pp = regionprops(mk, 'BoundingBox');

num_reg = length(pp);

if isempty(cc)
    cc = colormap(lines(num_reg));
end

use_lighting = 1;

sz = size(mk);

for ii = 1:num_reg
    
    [r,c,z] = getBoundingBox3(pp, ii, 11, sz);
    
    mkk = ~~((mk(r,c,z)==ii) + ~~surfs(r,c,z));
    mkk = imfilter3(double(mkk), [1,1,1]);
    
    mk_ = 0*mk;
    mk_(r,c,z) = mkk;
    
    cdata = cc(ii,:);
    
 fv = isosurface(double(mk_),3/4);
 p = patch(fv);

 set(p,...
    'facecolor', cdata,...
    'edgealpha', 0,...
    'facealpha', 1);
        
end

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
    
end

