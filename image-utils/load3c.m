function im =  load3(dd, dtype)

slash = filesep;
if ~exist('dtype', 'var')
    dtype = 'uint16';
end


if isdir(dd)
    
    ff = dir([dd slash '*.tif']);
    
    if isempty(ff)
        ff = dir([dd slash '*.tiff']);
    end
    
    switch dtype
        
        case 'double'
            for ii = 1:length(ff)
                
                im(:,:,ii,:) = double( imread([dd slash ff(ii).name]));
                
            end
        case 'uint8'
            for ii = 1:length(ff)
                
                im(:,:,ii,:) = uint8( imread([dd slash ff(ii).name]));
                
            end
        case 'uint16'
            for ii = 1:length(ff)
                
                im(:,:,ii,:) = uint16( imread([dd slash ff(ii).name]));
                
            end
    end
    
    
else
    
    numFrames = length(imfinfo(dd));
    
    switch dtype
        
        case 'double'
            for ii = 1:numFrames
                
                im(:,:,ii,:) = double( imread(dd, 'index', ii));
                
            end
        case 'uint8'
            for ii = 1:numFrames
                
                im(:,:,ii,:) = uint8( imread(dd, 'index', ii));
                
            end
        case 'uint16'
            for ii = 1:numFrames
                
                im(:,:,ii,:) = uint16( imread(dd, 'index', ii));
                
            end
    end
    
    
    
end

