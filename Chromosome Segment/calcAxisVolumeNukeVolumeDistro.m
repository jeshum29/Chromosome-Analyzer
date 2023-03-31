function [gauss_longest, gauss_shortest] = calcAxisLengthDistros(data)

for ii = 1:length(data)
    
    len = data(ii).props(12,:);
    vol = data(ii).props(13,:);
    
    nuke_vol(ii) = data(ii).props(9,1);
    
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

max_axis_volume = max(vols,[],2);
min_axis_volume = min(vols,[],2);
nuke_volume = nuke_vol';

C = calc_cov_matrix([max_axis_volume, nuke_volume]);
M = [mean(max_axis_volume), mean(nuke_volume)];

gauss_longest.M = M;
gauss_longest.C = inv(C);

C = calc_cov_matrix([min_axis_volume, nuke_volume]);
M = [mean(min_axis_volume), mean(nuke_volume)];

gauss_shortest.M = M;
gauss_shortest.C = inv(C);



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

