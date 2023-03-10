function [xxf,yyf,zzf,intensities] = trace3(im, xx, yy, zz, crop)

if length(crop)==1
    crop = [1,1,1]*crop;
end

for jj = 1:length(xx)-1
    
    positionVec(jj,:) = [xx(jj),yy(jj),zz(jj)];
    tangentVec(jj,:) = [xx(jj+1) - xx(jj), yy(jj+1) - yy(jj), zz(jj+1) - zz(jj)];
    
end

% use the tangentVec at the end-1 position for the last point
% this ensures xxf is the same length as xx

positionVec(length(xx),:) = [xx(end),yy(end),zz(end)];
tangentVec(length(xx),:) = tangentVec(length(xx)-1,:);


for jj = 1:size(positionVec,1)
    
    [immR, imm, tangentVecR] = getRotatedCenteredImage(im, im*0+1, crop, positionVec(jj,:), tangentVec(jj,:));
    
    image_center = (size(immR)+1)/2;
    immR_crop = immR(image_center(1) - 1:image_center(1) + 1, :, :);
    
    [X,Y,Z] = ndgrid(1:size(immR_crop,1), 1:size(immR_crop,2), 1:size(immR_crop,3));
    
    Y_CEN = sum(Y(:).*immR_crop(:))/sum(immR_crop(:));
    Z_CEN = sum(Z(:).*immR_crop(:))/sum(immR_crop(:));
    
    center_final = [image_center(1), Y_CEN, Z_CEN] - image_center;
        
%     if tangentVec(jj,2)==0
%         tangentVec(jj,2) = 1e-3;
%     end
        
    center_final = rotateVectorBack(center_final, tangentVec(jj,:));
    
    xxf(jj) = xx(jj) + center_final(1);
    yyf(jj) = yy(jj) + center_final(2);
    zzf(jj) = zz(jj) + center_final(3);
    
    intensities(jj) = immR(image_center(1), round(Y_CEN), round(Z_CEN));
        
    if jj>1
    
        dr = [xxf(end) - xxf(end-1), yyf(end) - yyf(end-1), zzf(end) - zzf(end-1)];
        step_size = sqrt(sum(dr.^2));
    
    
        if 0

            immR = autogain(immR);
            
            vizIm(immR - mean(immR(:)),2);
            vizIm(immR_last - mean(immR_last(:)), 22);

            '';
        end
        
        
        if 0
            
            figure(9);
            plot3(xxf,yyf,zzf);
            axis vis3d;
            
            if step_size > .7
                
            '';
            end
            
        end
        
    end
    

    immR_last = immR;
    immR_crop_last = immR_crop;
    
    
end



'';






