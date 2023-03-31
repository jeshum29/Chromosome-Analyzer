function batchWriteImages(rootDir, Z_RANGE, FILTER_STRING)

if ~exist('Z_RANGE', 'var')
    Z_RANGE = [];
end


% reconstructed OMX image pixel size
PX_SZ = [0.0649/2 0.0649/2 0.2]*1000; %in nm
RESZ_RATIO = PX_SZ(3)/PX_SZ(1);

slash = filesep;
fnameList = dir([rootDir slash '*' FILTER_STRING '*.tif']);

for ii = 1:length(fnameList)
    fnameList(ii).name = [rootDir filesep fnameList(ii).name];
end

SEARCH_STRING = [FILTER_STRING '_C'];

% generate a list of unique file names with the channel tag removed
counter = 0;
for ii = 1:length(fnameList)
    
    [path, fname, ext] = fileparts(fnameList(ii).name);
    fname_root = fname(1:max(length(SEARCH_STRING) - 1 + regexp(fname, SEARCH_STRING)));
        
    unique_fname = 1;
    
    if isempty(fname_root)
        continue;
    end
    
    for jj = ii+1:length(fnameList)
        
        [path, fname, ext] = fileparts(fnameList(jj).name);
        fname_root_ = fname(1:max(length(SEARCH_STRING) - 1 + regexp(fname, SEARCH_STRING)));
        
        if strcmpi(fname_root, fname_root_) || isempty(fname_root)
            unique_fname = 0;
        end
    end
    
    if unique_fname
        counter = counter + 1;
        fnameListUnique(counter).name = fname_root;
    end
    
end

try
    fnameListUnique;
catch
    ''
end

if ~exist('fnameListUnique', 'var')
    ['WARNING: no DVs found in ' rootDir]
    return;
end

% write images for each file
for ii = 1:length(fnameListUnique)
    
    fname_root = fnameListUnique(ii).name;
    fname_root = fname_root(1:end-2);
    
    writeDir = [rootDir slash fname_root]
        
    chans = [];
    
    for jj = 1:length(fnameList)
        
        [path, fname, ext] = fileparts(fnameList(jj).name);
        
        fname_root_ = fname(1:max(regexp(fname, 'C')));
        
        if strcmpi(fnameListUnique(ii).name, fname_root_)
            chans = [chans, str2double(fname(max(regexp(fname,'C'))+1))];
        end
                
    end
    
    chans = sort(chans);
    for jj = 1:length(chans)
        
        if ~exist([rootDir slash fname_root '_C' num2str(chans(jj)) '.tif'])
            ''
        end
        
        if exist([writeDir slash 'im' num2str(chans(jj)) '.tif'])
            continue;
        end
        
        im = load3([rootDir slash fname_root '_C' num2str(chans(jj)) '.tif']);
        sz = size(im);
        
        % if the TIFF is actually a projection, we should skip
        if length(sz)==2
            continue;
        end
        
        % if we've gotten this far, make the subdirectory and write the
        % resampled images 
        if ~isdir(writeDir)
            mkdir(writeDir);
        end
        
%         % crop in z if the stack is much too large (more than 7um)
%         MIN_SZ = 999;
%         
%         if jj==1 
%             if sz(3) > MIN_SZ && isempty(Z_RANGE)
%                 
%                 prof       = calcIntensityProfileZ(im);
%                 [pp, ind]  = max(prof);
%                 Z_RANGE    = max(1,ind-20):min(size(im,3),ind+20);
%             else
%                 Z_RANGE = 1:sz(3);
%             end
%         end
    
        Z_RANGE = 1:sz(3);
        
        im = im(:,:, Z_RANGE);
        
        sz = size(im);
        im = imresize3(im, [sz(1), sz(2), sz(3)*RESZ_RATIO]);

        save3(uint16(im), [writeDir slash ], ['im' num2str(chans(jj))], 0);
        
        % check that the full image was written (this compensates for 
        % an apparent bug in imwrite)
        im_ = load3([writeDir slash 'im' num2str(chans(jj)) '.tif']);
        if size(im_, 3) ~= size(im,3)
            save3(uint16(im), [writeDir slash ], ['im' num2str(chans(jj))], 0);
        end
        
        % write a projection
        imwrite(autogain(max(im,[],3)), [writeDir slash 'im' num2str(chans(jj)) '_PROJ.tif']);

    end
end

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

% r = autogain(r);
% g = autogain(g);
% b = autogain(b);

r = r(:,:,zr); g = g(:,:,zr); b = b(:,:,zr);

imR = squeeze(max(r_scale*(r - mean(r(:))*r_off), [], DIM));
imG = squeeze(max(g_scale*(g - mean(g(:))*g_off), [], DIM));
imB = squeeze(max(b_scale*(b - mean(b(:))*b_off), [], DIM));

imR  = autogain(double(imR).^(r_gamma));
imG  = autogain(double(imG).^(g_gamma));
imB  = autogain(double(imB).^(b_gamma));

im = cat(3,autogain(imR),autogain(imG),autogain(imB));


end

