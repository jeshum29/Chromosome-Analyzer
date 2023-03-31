function gauss_axis_vols = calcAxisLengthDistros(data)

for ii = 1:length(data)
    
    len = data(ii).props(12,:);
    vol = data(ii).props(13,:);
    
    nuke_vol = data(ii).props(11,1);
    
    vol = vol(~~vol);
    len = len(~~len);
    
    
    vols(ii,1:length(vol)) = vol;
    lens(ii,1:length(len)) = len;
    
    
end

% sort by size
vols = sort(vols,2,'descend');
lens = sort(lens,2,'descend');

% only look at the six longest
vols = vols(:, 1:6);
lens = lens(:, 1:6);

% normalize volumes by total volume in each nucleus
vols = vols ./ repmat(sum(vols,2), 1, 6);

% one volume is now redundant
vols = vols(:,1:5);

C = calc_cov_matrix(vols);
vol_mean = mean(vols,1);

gauss_axis_vols.M = vol_mean;
gauss_axis_vols.C = inv(C);

cc = loadColors;

% for ii = 1:5
%     
%     vol = vols(:,ii) - vols(:,ii+1);
%     vol = vol(~~vol);
%     
%     [h,c] = hist(vol,30);
%     
%     plot(c,h, 'color', cc(ii,:));
%     hold on
%     
% end

% c = 0;
% for ii = 1:size(vols,1)
%     for jj = ii:size(vols,1)
%         c = c+1;
%         dd(c) = sqrt(sum((vols(ii,:) - vols(jj,:)).^2));
%     end
% end
% 
% figure;
% [h,c] = hist(dd,100);
% plot(c,h);



% [vecs,vals] = eig(calc_cov_matrix(vols));
%
% vol_proj = (sum((vols - repmat(vol_mean, size(vols,1), 1)) .* repmat(vecs(:,end)', size(vols,1), 1), 2));
% [h,c] = hist(vol_proj, 30);
% plot(c,h);


% for ii = 1:size(vols,1)
%     vol = vols(ii,:);
%     prob(ii) = exp( -1/2 * (vol - vol_mean) * (C) * (vol - vol_mean)' );
% end



'';

end



function covMatrix = calc_cov_matrix(vol)

vol = vol - repmat(mean(vol,1), size(vol,1), 1);

for ii = 1:size(vol,2)
    for jj = 1:size(vol,2)
        covMatrix(ii,jj) = sum(vol(:,ii).*vol(:,jj));
    end
end

covMatrix = covMatrix/size(vol,1);


end

