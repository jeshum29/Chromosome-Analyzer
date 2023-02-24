function [gauss_longest, gauss_shortest] = calcAxisLengthDistros(data)

for ii = 1:length(data)
    
    len = data(ii).props(11,:);
    vol = data(ii).props(12,:);
    
    total_len(ii) = sum(len);
    
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

max_axis_len = max(lens,[],2);
min_axis_len = min(lens,[],2);
total_len = total_len';

C = calc_cov_matrix([max_axis_len, total_len]);
M = [mean(max_axis_len), mean(total_len)];

gauss_longest.M = M;
gauss_longest.C = inv(C);

C = calc_cov_matrix([min_axis_len, total_len]);
M = [mean(min_axis_len), mean(total_len)];

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

