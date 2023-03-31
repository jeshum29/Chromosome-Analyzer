
function [vecs, vals, cen] = calcPrincipleAxes(im)

sz = size(im);
im = double(im);

[X,Y,Z] = ndgrid(1:sz(1), 1:sz(2), 1:sz(3));


im_sum = sum(im(:));

Cx = sum(X(:).*im(:))/im_sum;
Cy = sum(Y(:).*im(:))/im_sum;
Cz = sum(Z(:).*im(:))/im_sum;

cen = [Cx, Cy, Cz];

X = X - Cx;
Y = Y - Cy;
Z = Z - Cz;

Ixy = sum(X(:) .* Y(:) .* im(:));
Ixz = sum(X(:) .* Z(:) .* im(:));
Iyz = sum(Z(:) .* Y(:) .* im(:));

Ixx = sum(X(:) .* X(:) .* im(:));
Iyy = sum(Y(:) .* Y(:) .* im(:));
Izz = sum(Z(:) .* Z(:) .* im(:));

I = [...
    [Ixx, Ixy, Ixz];...
    [Ixy, Iyy, Iyz];...
    [Ixz, Iyz, Izz];...
    ];

I = I/im_sum;


[vecs, vals] = eig(I);

vals = diag(vals);

end
