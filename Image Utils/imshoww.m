function imshoww(im, minmax)

if ~exist('minmax', 'var')
    minmax = [];
end

imshow(im,[minmax],'initialmagnification', 'fit');

end
