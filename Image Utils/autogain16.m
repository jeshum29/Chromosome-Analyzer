function im = autogain( im)

im = double(im);

im_min = min(im(:));
im = im - im_min;

im_max = max(im(:));
im = uint16(65535*im/im_max);