function [nd] = spotAxis3(nd, SPOT_CHAN, THRESH, redo, SMOOTHING)

if ~exist('redo')
    redo = 0;
end

if ~exist('SMOOTHING')
    SMOOTHING = 3;
end

if ~isfield(nd, 'PX_SZ')
    PX_SZ = 110;
    'Warning: no pixel size specified'
else
    PX_SZ = nd.PX_SZ;
end


eval(['im = nd.im' num2str(SPOT_CHAN) ';']);

if PX_SZ > 100

    im = double(imresize3(im, size(im)*2));
    mk = double(imresize3(imdilate(nd.mk_back, ones(3,3,3)), size(nd.mk_back)*2));

else
    
    im = double(im);
    mk = double(imdilate(nd.mk_back, ones(3,3,3)));
    
end

% --------------------------------------------------------------
%
% removal of mask mk from nd.mk_back
% this mask has, as far as I know, always been used, and
% has certainly been used for all RAD-51 OMX data.
% 
% It has the effect of sometimes masking real RAD-51 foci,
% so it should not be used.
%
% The following line added (i.e., use of mk_back removed)
% on 2015_02_13  
% -KCC
%

mk = 0*mk + 1;

%
% --------------------------------------------------------------


if redo
    
    spots = spotter3(im, mk, THRESH, SMOOTHING);
    eval(['nd.spots_im' num2str(SPOT_CHAN) ' = spots;']);
    
else
    spots = nd.spots;
end

if isempty(spots)
    'warning! no spots found!'
    return;
end

for ii = 1:length(spots)
    
    score(ii) = spots(ii).intensity_score;
    
end

% sd = nd.sdata;
% 
% [s,ord] = sort(score, 'descend');
% 
% for ii = 1:length(spots)
%     
%     sp = spots(ord(ii));
%     
%     for jj = 1:length(sd)
%         
%         if ~isempty(sd(jj).trace)
%             spot_pos = repmat(sp.r, size(sd(jj).trace,1), 1);
%             dd = sqrt(sum((sd(jj).trace - spot_pos).^2,2));
%             delta(jj) = min(dd);
%         else
%             delta(jj) = 2*MAX_AXIS_SPOT_DISTANCE;
%         end
%     end
%     
%     [d, jj_min] = min(delta);
%     
%     if d < MAX_AXIS_SPOT_DISTANCE
%         nd.spots(ord(ii)).axis_id = jj_min;
%     else
%         nd.spots(ord(ii)).axis_id = 0;
%     end
%     
%     
% end
% 
% 
% for ii = 1:length(nd.sdata)
%     
%     spot_ids = [];
%     
%     for jj = 1:length(nd.spots)
%         if nd.spots(jj).axis_id==ii
%             spot_ids = [spot_ids, jj];
%         end
%     end
%     
%     if isempty(spot_ids)
%         
%     else
%         
%         ints = [nd.spots(spot_ids).intensity_score];
%         [m, ind] = max(ints);
%         
%         nd.sdata(ii).spot_id = spot_ids(ind);
%         
%         for jj = 1:length(spot_ids)
%             
%             if spot_ids(jj)~=spot_ids(ind)
%                 nd.spots(spot_ids(jj)).axis_id = 0;
%             end
%             
%         end
%         
%     end
%     
% end
% 


%
%
% vizIm(cat(4,...
%     autogain(nd.im1 - mean(nd.im1(:))), ...
%     autogain(nd.im2 - mean(nd.im2(:))),...
%     0*nd.im1), 2);


