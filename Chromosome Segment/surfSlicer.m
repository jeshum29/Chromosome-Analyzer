function surfSlicerData = surfSlicer(nd)

DIST = 6;
sz = size(nd.surfL);
crop = [1,1,1]*8 + DIST + 1;

surfSlicerData.DIST = DIST;
surfSlicerData.step_size = 1;
surfSlicerData.crop = crop;


for ii = 1:length(nd.ppS)
    
    clear tangentVec positionVec reg_kymo surf_kymo
    
    regs = nd.sd{ii}.neighbors;
    
    pos = round(nd.ppS(ii).Centroid);
    
    D = DIST*2;
    
    rr = (-D:D) + pos(2);
    cc = (-D:D) + pos(1);
    zz = (-D:D) + pos(3);
    
%     [rr,cc,zz] = getBoundingBox3(nd.ppR, regs, 3, sz);
%     [rrr,ccc,zzz] = getBoundingBox3(nd.ppS, ii, 3, sz);
%     
%     rr = min([rr,rrr]):max([rr,rrr]);
%     cc = min([cc,ccc]):max([cc,ccc]);
%     zz = min([zz,zzz]):max([zz,zzz]);
    
    rr = max(rr(1), 1):min(rr(end), sz(1));
    cc = max(cc(1), 1):min(cc(end), sz(2));
    zz = max(zz(1), 1):min(zz(end), sz(3));
    
    % neighborRegs = (nd.regL(rr,cc,zz)==regs(1)) + (nd.regL(rr,cc,zz)==regs(2)) + (nd.surfL(rr,cc,zz)==ii);
        
    surff = nd.surfL(rr,cc,zz)==ii;
    regg = ~~(nd.regL(rr,cc,zz) + nd.surfL(rr,cc,zz));
    
    reggL = bwlabeln(regg, 6);
    ind = mode(reggL(surff(:)));

    for jj = 1:max(reggL(:))
        if ~(jj==ind)
            regg(reggL==jj) = 0;
        end
    end
    

    regg_ = zeros(size(regg)+crop*2); 
    regg_(crop(1)+1:end-crop(1), crop(2)+1:end-crop(2), crop(3)+1:end-crop(3)) = regg;
    regg = regg_; 
    
    surff_ = zeros(size(surff)+crop*2); 
    surff_(crop(1)+1:end-crop(1), crop(2)+1:end-crop(2), crop(3)+1:end-crop(3)) = surff;
    surff = surff_; 
    
    [vecs, vals, cen] = calcPrincipleAxes(surff);
    
    [vv, min_ind] = min(vals);
    ortho_vec = vecs(:, min_ind)';
    
    steps = -surfSlicerData.DIST:surfSlicerData.step_size:surfSlicerData.DIST;
    
    for jj = 1:length(steps)
        
        positionVec(jj,:) = cen + ortho_vec*steps(jj);
        
    end
    
    tangentVec = positionVec(2:end,:) - positionVec(1:end-1,:);
    positionVec = positionVec(2:end,:);
    
    mask = ones(size(regg));
        
    for jj = 1:size(positionVec,1)

        [immR, imm] = getRotatedCenteredImage(regg, mask, [10,10,10], positionVec(jj,:), tangentVec(jj,:));
        image_center = (size(immR)+1)/2;
        reg_kymo(jj,:,:) = immR(image_center(1),:,:);
        
        [immR, imm] = getRotatedCenteredImage(surff, mask, [10,10,10], positionVec(jj,:), tangentVec(jj,:));
        image_center = (size(immR)+1)/2;
        surf_kymo(jj,:,:) = immR(image_center(1),:,:);
        
    end
    
    surfSlicerData.reg_kymo(ii).kymo = reg_kymo;
    surfSlicerData.surf_kymo(ii).kymo = surf_kymo;
    
    
    
    if 0
        
        vizVolume3(regg, [1,0,0],1,1,1,1/2);
        vizVolume3(surff, [0,1,0], 1, 0, 0, 1/2);

        reg_sum = sum(sum(reg_kymo,3),2);
        surf_sum = sum(sum(surf_kymo,3),2);

        vizIm(reg_kymo,2);
        figure(22);clf;
        hold on
        plot([sum(sum(reg_kymo,3),2),sum(sum(surf_kymo,3),2)]);


        [mx, surf_pos] = max(surf_sum);

        reg_sum = reg_sum(max(1,surf_pos-4):min(end,surf_pos+4));
        reg_dd = (reg_sum(2:end) - reg_sum(1:end-1))./reg_sum(2:end);

        reg_dd(reg_sum(2:end)==0) = 0;

        midd = round(length(reg_dd)/2);
        param = sum(reg_dd(midd:end)) - sum(reg_dd(1:midd))
    
    end
    
    '';
end

end


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
