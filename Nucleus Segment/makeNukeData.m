function N = makeNukeData(gonad, CROP, OVERWRITE_FLAG)


%-------------------------------------------------------------------------------------
%
% CROP SIZE IS IMPORTANT
%
%
%-------------------------------------------------------------------------------------

N = 0;

slash = filesep;
dd = gonad.writeDir;

MIN_NUKE_VOLUME  = gonad.MIN_NUKE_VOLUME;
MAX_NUKE_VOLUME  = gonad.MAX_NUKE_VOLUME;
MAX_NUKE_DIA     = gonad.MAX_NUKE_DIA;
MIN_NUKE_DIA     = gonad.MIN_NUKE_DIA;
MAX_BORDER_AREA  = gonad.MAX_BORDER_AREA;

try
  mkF = load3([dd slash 'mkF3_reseg.tif']);
catch
    'FAILURE: mkF3_reseg.tif does not exist'
    return;
end

for ii = 1:gonad.NUM_CHAN
    
    ii_str = num2str(ii);
    im{ii} = load3([dd slash 'im' ii_str '.tif']);
    
end

sz = size(im{1});
    
imb = imread([dd slash 'mkF3_reseg_color_proj.tif']);
    
props = regionprops(mkF, 'Area', 'BoundingBox');

if 0
    
    imR = imb(:,:,1);
    imG = imb(:,:,2);
    imB = imb(:,:,3);
    
    figure(66);
    imshow(cat(3, imR, imG, imB),[]);
    hold on
    
    FONT_SZ = round(48 * size(imR,1) / 1000);
    
    sz = size(mkF);
    
    props = regionprops(mkF, 'Area', 'BoundingBox','Centroid');
    
    for ii = 1:length(props)
        
        [rr,cc,zz] = getBoundingBox3(props, ii, 0, sz);
        dia = max([rr(end) - rr(1), cc(end) - cc(1), zz(end) - zz(1)]);
        
        include_flag = ...
            props(ii).Area > MIN_NUKE_VOLUME && ...
            props(ii).Area < MAX_NUKE_VOLUME && ...
            dia < MAX_NUKE_DIA;
        
        ss = num2str(ii);
        
        if include_flag
            
            text(...
                props(ii).Centroid(1),...
                props(ii).Centroid(2),...
                ss, 'fontsize', FONT_SZ,'fontweight', 'bold', 'color', [0,.8,0]);
            
        else
            
            text(...
                props(ii).Centroid(1),...
                props(ii).Centroid(2),...
                ss, 'fontsize', FONT_SZ,'fontweight', 'bold', 'color', [0,0,1]);
            
        end
    end
    
end


ndata = [];

for ii = 1:length(props)
    
    [rr,cc,zz] = getBoundingBox3(props, ii, 0, sz);
    dia = max([rr(end) - rr(1), cc(end) - cc(1), zz(end) - zz(1)]);
    
    [rr,cc,zz] = getBoundingBox3(props, ii, CROP, sz);
    
    rr = max(1, min(rr)) : min(sz(1), max(rr));
    cc = max(1, min(cc)) : min(sz(2), max(cc));
    zz = max(1, min(zz)) : min(sz(3), max(zz));
    
try
    mkk = mkF(rr,cc,zz)==ii;
catch
    ''
end
    
    include_flag = isNuke2(...
        props(ii).Area, dia, mkk,...
        MIN_NUKE_VOLUME,...
        MAX_NUKE_VOLUME,...
        MAX_NUKE_DIA,...
        MAX_BORDER_AREA);
    
    if isfield(gonad, 'deleted_regions')
        if ~isempty(intersect(gonad.deleted_regions, ii))
            include_flag = 0;
        end
    end
    
    if include_flag
        
        try
            mkkF = mkF(rr,cc,zz);
            add_nuke = 1;
        catch
            add_nuke = 0;
            'region too close to edge - skipping'
        end
        
        if add_nuke
            ndata{end+1}.label = ii;
            ndata{end}.rr = rr;
            ndata{end}.cc = cc;
            ndata{end}.zz = zz;
            ndata{end}.props = props(ii);
        end
        
    end
    
end

if OVERWRITE_FLAG && isdir([dd slash 'nukes'])

    'WARNING in makeNukeData: Overwriting nucleus data!'
    
elseif ~isdir([dd slash 'nukes'])
   
    % safe to write
    
else
    
    'WARNING in makeNukeData: no overwrite permission; exiting.'
    return;
    
end
    


try
    rmdir([dd slash 'nukes'],'s');
end

try
    rmdir([dd slash 'gallery'],'s');
end

if isfield(gonad, 'nuke_ids')
    gonad = rmfield(gonad, 'nuke_ids');
end
if isfield(gonad, 'finalized_nukes')
    gonad = rmfield(gonad, 'finalized_nukes');
end

save( [gonad.writeDir filesep 'gonad.mat'], 'gonad');

eval(['!del ' dd slash 'ndata.mat']);
save([dd slash 'ndata.mat'], 'ndata');

ndata = load([dd slash 'ndata.mat']);
ndata = ndata.ndata;

mkdir([dd slash 'gallery']);
mkdir([dd slash 'nukes']);

N = length(ndata);



for ii = 1:length(ndata)
    
    nd = ndata{ii};
    
    fn = [num2str(nd.label), '_'];
    
    ddd = [dd slash 'nukes' slash num2str(nd.label) slash];
    
    mkk = mkF(nd.rr,nd.cc,nd.zz)==nd.label;
    mkk = imfill(mkk, 'holes');
    
    save3(autogain(mkk), ddd, [fn 'mkk'], 0);
    writeTIFFmeta([ddd slash fn 'mkk.tif'], 0, 1, 1);
    
    nd.mk = logical(mkk);
    nd.mk_back = imdilate(mkk, ones(3,3,3));
    
    
    for jj = 1:gonad.NUM_CHAN
        
        imm = im{jj};

        imm = imm(nd.rr,nd.cc,nd.zz);

        eval(['nd.im' num2str(jj) '= imm;']);
        
        save3(autogain(imm), ddd, [ 'im' num2str(jj)], 0);
        writeTIFFmeta([ddd slash 'im' num2str(jj) '.tif'], 0, max(imm(:)), 1);
        
        imm = imm; %.*uint16(nd.mk_back);
        imPRJ = makeProjection(imm,imm,imm, 3);
        
        imPRJ = imresize(imPRJ, 6, 'nearest');
        
        imwrite(autogain(imPRJ), [ddd slash 'im' num2str(jj) '_PRJ.tif']);
        imwrite(autogain(imPRJ), [dd slash 'gallery' slash num2str(ii) '_' num2str(jj) '_PRJ.tif']);
        
    end
    
    nd.dir = ddd;
    save([ddd 'nd'], 'nd');
    
    
end


function im = makeProjection(r,g,b, DIM)

zr = 1:ceil(size(r,3)/1);

sc = 2;

r_scale = sc;
g_scale = sc;
b_scale = sc;

r_off = 1;
g_off = 1;
b_off = 1;

ga = .8;

r_gamma = ga;
g_gamma = ga;
b_gamma = ga;

r = r(:,:,zr); g = g(:,:,zr); b = b(:,:,zr);


imR = squeeze(max(r_scale*(r - mean(r(:))*r_off), [], DIM));
imG = squeeze(max(g_scale*(g - mean(g(:))*g_off), [], DIM));
imB = squeeze(max(b_scale*(b - mean(b(:))*b_off), [], DIM));

imR  = autogain(double(imR).^(r_gamma));
imG  = autogain(double(imG).^(g_gamma));
imB  = autogain(double(imB).^(b_gamma));


im = cat(3,autogain(imR),autogain(imG),autogain(imB));


