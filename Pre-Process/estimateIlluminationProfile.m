function estimateIlluminationProfile(dds)

BACKGROUND = 5000;


for ii = 1:length(dds)
    
    
    im_ = imread([dds(ii).name filesep 'im1_PROJ.tif']);
    
    im_(im_<BACKGROUND) = 0;
    
    im(:,:,ii) = im_;
    
    imF(:,:,ii) = imfilter(im(:,:,ii), fspecial('gaussian', 45, 15), 'same');
    
end


angles = 0:5:90;

cen = size(im,1)/2;
crop = 10;

for ii = 1:size(im,3)
    for jj = 1:length(angles)
        
        
        imR = imrotate(imF(:,:,ii), angles(jj),'bilinear','crop');
        
        sliver = imR(cen-crop:cen+crop,:);        
        profH(ii,jj,:) = sum(sliver,1) ./ (sum(~~sliver,1)+1);
        profH(ii,jj,:) = max(sliver,[],1);
        
        sliver = imR(:, cen-crop:cen+crop)';
        profV(ii,jj,:) = sum(sliver,1) ./ (sum(~~sliver,1)+1);
        profV(ii,jj,:) = max(sliver,[],1);
        
    end
end


prof = cat(2, profH, profV);

profMean = sum(sum(prof,1),2) ./ sum(sum(~~prof,1),2);


plot(squeeze(mean(mean(cat(2, profH, profV),1),2)))

'';