function imF = imfilter3(im, gaussian_sz)


sz = size(im);

pad = max(1,2*max(ceil(gaussian_sz)));

im_ = zeros(sz(1) + 2*pad, sz(2) + 2*pad, sz(3));

im_(pad + 1 : sz(1) + pad , ...
    pad + 1 : sz(2) + pad , ...
    1 : sz(3))...
    = im;

im_(1:pad, :, :) = repmat(im_(pad + 1,:,:), [pad, 1,1]);
im_(pad + sz(1) + 1 : end, :, :) = repmat(im_(sz(1) - pad,:,:), [pad, 1, 1]);


im_(:, 1:pad, :) = repmat(im_(:, pad + 1,:), [1, pad, 1]);
im_(:, pad + sz(2) + 1 : end, :) = repmat(im_(:, sz(2) - pad, :), [1, pad, 1]);


 kernel = gaussian3(size(im_), gaussian_sz);

 imF = ifftshift(ifftn(fftn(im_).*fftn(kernel)));


imF = imF(pad+1:sz(1)+pad, pad+1:sz(2)+pad, 1:sz(3));

