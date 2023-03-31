function nd = makeFocusData(nd)

% list of image channels corresponding to the two possible focus channels
focusChannel =     nd.focusChannel;

% focus-axis identities
focusAxisTable =   nd.focusAxisTable;

if focusChannel(1)
    try
        eval(['spots_1 = nd.spots_im' num2str(focusChannel(1)) ';']);
        spots_1_chan = focusChannel(1);
    catch
        ['no spots in channel ' num2str(focusChannel(1))]
        spots_1 = [];
    end
else
    spots_1 = [];
end

if focusChannel(2)
    try
        eval(['spots_2 = nd.spots_im' num2str(focusChannel(2)) ';']);
        spots_2_chan = focusChannel(2);
    catch
        ['no spots in channel ' num2str(focusChannel(2))]
        spots_2 = [];
    end
else
    spots_2 = [];
end

% ============================================================--
%
%  here, we list each focus and assume that each corresponds to a different
%  chromosome. 
%
% ============================================================--

focusPositions =  [];
focusID =         [];
focusChannel =    [];

cc = 0;

if ~isempty(spots_1)
    
    clear score
    for ii = 1:length(spots_1)
        score(ii) = spots_1(ii).intensity_score;
    end

    [ss, ord] = sort(score, 'descend');
    chrID =     focusAxisTable(:, 1);
    numSpots =  length(spots_1);
    
    ind = find(~~chrID);
    
    for ii = 1:min(length(ind), numSpots)
        cc = cc+1;
        fd{cc}.position =  spots_1(ord(ind(ii))).r;
        fd{cc}.chrID =     chrID(ind(ii));
        fd{cc}.channel =   spots_1_chan;
        fd{cc}.index =     ord(ind(ii));
    end
    
end

if ~isempty(spots_2)
    
    clear score
    for ii = 1:length(spots_2)
        score(ii) = spots_2(ii).intensity_score;
    end

    [ss, ord] = sort(score, 'descend');
    chrID =     focusAxisTable(:, 2);
    numSpots =  length(spots_2);
    
    ind = find(~~chrID);
    
    for ii = 1:min(length(ind), numSpots)
        cc = cc+1;
        fd{cc}.position =  spots_2(ord(ind(ii))).r;
        fd{cc}.chrID =     chrID(ind(ii));
        fd{cc}.channel =   spots_2_chan;
        fd{cc}.index =     ord(ind(ii));
    end
    
end

nd.fd = fd;


