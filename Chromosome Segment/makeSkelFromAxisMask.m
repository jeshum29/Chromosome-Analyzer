function sk = makeSkelFromAxisMask(mk,spur_len)

mk = ~~mk;

% this should remove holes
mk = imfill(mk,'holes');

pp = regionprops(bwlabeln(~mk, 6),'PixelIdxList');

% just in case, remove holes again
for ii = 1:length(pp)
    if length(pp(ii).PixelIdxList) < 1000
        mk(pp(ii).PixelIdxList) = 1;
    end
end

[mkL, num] = bwlabeln(mk, 6);

pp = regionprops(mkL, 'BoundingBox');
sz = size(mk);

sk = 0*mk;

for ii = 1:num
    
    [rr,cc,zz] = getBoundingBox3(pp, ii, 1, sz);
    
    mkk = mkL(rr,cc,zz)==ii;
    skk = ~~skel3(double(mkk));
    
    sk_ = sk(rr,cc,zz);
    sk_(skk) = ii;
    sk(rr,cc,zz) = sk_;
    
end

spurs = (~~sk) - deleteSpurs(~~sk, spur_len);

sk(~~spurs) = 0;

