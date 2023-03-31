function writeTIFFfromDV(dvfile, writeProjection, writeStack)

meta = imreadBFmeta(dvfile);

% meta.width : image width
% meta.height : image height
% meta.zsize : number of z slices
% meta.nframes : number of time frames
% meta.channels : number of channels

NUM_Z = meta.zsize;
NUM_CHAN = meta.channels;

[dname, fname] = fileparts(dvfile);


for chan = 1:NUM_CHAN
    
    im = imreadBF(dvfile, 1:NUM_Z, 1, chan);

    if writeProjection
        if exist([dname filesep fname '_PROJ_C' num2str(chan) '.tif'])
            continue;
        end
        
%         im_ = im;
%         im_(im<0) = im(im<0) + 2^16 - 1;
%         im = im_;
            
        im_proj = max(im,[],3);
        im_proj = autogain(im_proj);
        imwrite(im_proj, [dname filesep fname '_PROJ_C' num2str(chan) '.tif'], 'compression', 'none');
    
    end
    
    if writeStack
        if exist([dname filesep fname '_C' num2str(chan) '.tif'])
            % continue;
        end
        im = autogain16(im);
        sz = size(im,3);
        
        save3(im, dname, [fname '_C' num2str(chan)], 0);
        
        im_ = load3([dname filesep fname '_C' num2str(chan) '.tif']);
        if size(im_,3)~=sz
            save3(im, dname, [fname '_C' num2str(chan)], 0);
        end
        
    end
end



end

function im = autogain16(im)

im = double(im);

im_min = min(im(:));
im = im - im_min;

im_max = max(im(:));
im = uint16(65535*im/im_max);

end

function im = autogain( im)

im = double(im);

im_min = min(im(:));
im = im - im_min;

im_max = max(im(:));
im = uint8(255*im/im_max);
end
