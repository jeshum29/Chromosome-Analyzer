function nd = vizAxisTrace(nd, figNum, clfFlag, cc, SPOT_CHAN)

try
    
   eval(['spots = nd.spots_im' num2str(SPOT_CHAN) ';']);
   
catch
    
    spots = nd.spots;
    
end


if ~exist('figNum', 'var')
    figNum = 99;
end

if ~exist('clfFlag','var')
    clfFlag = 1;
end

if ~exist('cc','var')
    cc = ones(length(spots), 3);
end

if ~isfield(nd, 'PX_SZ')
    PX_SZ = 110;
    'Warning: no pixel size specified'
else
    PX_SZ = nd.PX_SZ;
end


if PX_SZ > 100
    scale = 1/2;
else
    scale = 1;
end

spot_offset = -1/2;

if size(cc, 1)==1
    cc = repmat(cc, length(spots), 1);
end

ff = figure(figNum);

if clfFlag
    clf;
end

set(ff, 'color', [0,0,0]);
nuke_cen = mean(reshape([spots.r],3,length([spots.r])/3)', 1);


for jj = 1:length(spots)
    
    hold on
    
    if isfield(spots(jj), 'axis_id') && spots(jj).axis_id
        ccc = cc(nd.sdata(spots(jj).axis_id).label,:);
    else
        ccc = cc(1,:);
    end
    
    scatter3(...
        spots(jj).r(2) * scale + spot_offset,...
        spots(jj).r(1) * scale + spot_offset,...
        spots(jj).r(3) * scale + spot_offset,...
        200, ccc, 'o', 'filled');
    
    text_offset = 5;
    unit_vec = spots(jj).r - nuke_cen;
    unit_vec = unit_vec ./ sqrt(sum(unit_vec.*unit_vec));
    
    if sum(isnan(unit_vec))
        unit_vec = [1,1,1]/sqrt(3);
    end
    
    text_pos = spots(jj).r*scale + text_offset * unit_vec;
    
    text_str = [...
        num2str(round(spots(jj).intensity_score*100)/100), ', ',...
        ];
    
    text(...
        text_pos(2), ...
        text_pos(1), ...
        text_pos(3), ...
        text_str, 'fontsize', 12, 'color', ccc);
    
end


daspect([1 1 1]);
view(3);
axis tight;
axis vis3d;
lighting none;

rotate3d;

set(gca, 'Position', [0,0,1,1]);

axis off



