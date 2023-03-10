function mkOut = makeColorMask(mk)

% mk = bwlabeln(mk, 26);

cc = autogain(colormap(parula(max(mk(:)))));

cc = autogain(colormap(hsv(max(mk(:)))));

[rr, ord] = sort(rand(1, size(cc,1)));

cc(:,1) = cc(ord,1);
cc(:,2) = cc(ord,2);
cc(:,3) = cc(ord,3);

mkR = uint8(0*mk); mkG = uint8(0*mk); mkB = uint8(0*mk);

props = regionprops(mk, 'Area','BoundingBox');

sz = size(mk);

for ii = 1:max(mk(:))
    
    if props(ii).Area < 500
        continue;
    end
    
    [r,c,z] = getBoundingBox3(props, ii, 3, sz);
    
    mkk = mk(r,c,z)==ii;
    
    
    mkRR = mkR(r,c,z);
    mkGG = mkG(r,c,z);
    mkBB = mkB(r,c,z);
    
    mkRR(mkk) = cc(ii,1);
    mkGG(mkk) = cc(ii,2);
    mkBB(mkk) = cc(ii,3);
    
    mkR(r,c,z) = mkRR;
    mkG(r,c,z) = mkGG;
    mkB(r,c,z) = mkBB;
   
    
end

mkOut = cat(4, mkR, mkG, mkB);