function imageEnergy = calcImageEnergy(im, opts)

Ix = imageDerivatives3D(im,opts.energySigma,'x');
Iy = imageDerivatives3D(im,opts.energySigma,'y');
Iz = imageDerivatives3D(im,opts.energySigma,'z');

im_line =  imfilter3(im, opts.energySigma);

im_edge = sqrt(Ix.^2 + Iy.^2 + Iz.^2);

imageEnergy = (-opts.edgeWeight * im_edge + opts.lineWeight * im_line);


end