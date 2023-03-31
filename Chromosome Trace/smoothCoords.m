

function [xx,yy,zz] = smoothCoords(x,y,z,interpsPerPoint,w)


ds = 1/interpsPerPoint;

tt = [x;y;z];

kk = fspecial('gaussian', 5*w*interpsPerPoint,w*interpsPerPoint);
kk_sz = round(size(kk,1)/2);
kk = kk(:,kk_sz);
kk = kk/sum(kk);

LL = size(tt,2);
LL_ = length(1:ds:LL);

for ii = 1:3

    x = tt(ii,:);
    x = interp1(1:LL, x, 1:ds:LL);

    x_pad = zeros(1,LL_+2*kk_sz) + x(1);
    x_pad(kk_sz:LL_+kk_sz-1) = x;
    x_pad(LL_+kk_sz:end) = x(end);

    x = conv2(x_pad, kk', 'same');
    x = x(kk_sz:LL_+kk_sz-1);

    if ii ==1
        xx = x;
    elseif ii==2
        yy = x;
    else
        zz = x;
    end

end
