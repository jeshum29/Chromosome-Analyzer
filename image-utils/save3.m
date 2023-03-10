function save3(im, dd, froot, saveEachFrame)

% function save3(im, dd, froot, saveEachFrame)


if ~exist('froot','var')
    froot = 'image';
end

if ~isdir(dd)
    mkdir(dd);
end

num_digits = length(num2str(size(im,3)));

fstr = ['%0' num2str(num_digits) 'd'];

for ii = 1:size(im,3)
    
    ff{ii} = [froot sprintf(fstr, ii) '.tif'];
    
end

slash = filesep;

if ~saveEachFrame
    
    eval(['!del ' dd froot '.tif']);
    
    for ii = 1:length(ff)
        try
        imwrite((squeeze(im(:,:,ii))), [dd slash froot '.tif'], 'compression', 'none', 'writemode', 'append');
        catch
            ''
        end
    end
    
else
    
    for ii = 1:length(ff)
        imwrite((squeeze(im(:,:,ii))), [dd slash ff{ii}], 'compression', 'none');
    end
    
end


