function im = autogain( im)

im = double(im);

im_min = min(im(:));
im = im - im_min;

im_max = max(im(:));
im = uint8(255*im/im_max);
