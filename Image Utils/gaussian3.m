function h = gaussian3(sz, std)

if length(std)==1
    std = [std,std,std];
elseif length(std)~=3
    'bad size vector'
    return;
end

if length(sz)==1
    sz = [sz,sz,sz];
elseif length(sz)~=3
    'bad size vector'
    return;
end

sz   = (sz-1)/2;
std = std.^2;

[x,y,z] = ndgrid(-sz(1):sz(1), -sz(2):sz(2), -sz(3):sz(3));
arg   = -(x.*x/std(1) + y.*y/std(2) + z.*z/std(3))/2;

h     = exp(arg);
h(h<eps*max(h(:))) = 0;

sumh = sum(h(:));

if sumh ~= 0,
    h  = h/sumh;
end