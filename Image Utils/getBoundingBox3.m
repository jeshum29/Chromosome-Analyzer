function [ii,jj, kk] = getBoundingBox(props, regions, pad, image_size)

if ~exist('pad', 'var')
    pad = 0;
end

xmin = 1;
ymin = 1;
zmin = 1;
ymax = image_size(1);
xmax = image_size(2);
zmax = image_size(3);

yyminmin = 9999;
xxminmin = 9999;
zzminmin = 9999;
yymaxmax = 0;
xxmaxmax = 0;
zzmaxmax = 0;

for nn = 1:length(regions)
    
    yymin = floor(props(regions(nn)).BoundingBox(2)) - pad;
    yymax = yymin + ceil(props(regions(nn)).BoundingBox(5)) + 1 + 2*pad;
    xxmin = floor(props(regions(nn)).BoundingBox(1)) - pad;
    xxmax = xxmin + ceil(props(regions(nn)).BoundingBox(4)) + 1 + 2*pad;
    zzmin = floor(props(regions(nn)).BoundingBox(3)) - pad;
    zzmax = zzmin + ceil(props(regions(nn)).BoundingBox(6)) + 1 + 2*pad;
    
    yymin = max([ymin,yymin]);
    xxmin = max([xmin,xxmin]);
    yymax = min([ymax,yymax]);
    xxmax = min([xmax,xxmax]);
    zzmin = max([zmin, zzmin]);
    zzmax = min([zmax, zzmax]);
    
    yyminmin = min([yyminmin,yymin]);
    xxminmin = min([xxminmin,xxmin]);
    yymaxmax = max([yymaxmax,yymax]);
    xxmaxmax = max([xxmaxmax,xxmax]);
    zzmaxmax = max([zzmaxmax, zzmax]);
    zzminmin = min([zzminmin, zzmin]);
    
end

ii = yyminmin:yymaxmax;
jj = xxminmin:xxmaxmax;
kk = zzminmin:zzmaxmax;

