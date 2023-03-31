function scores = batch(vols)


dataGlobal = loadGlobalData;
gauss_axis_vols = calcAxisVolumeDistro(dataGlobal);



for ii = 1:size(vols,1)
    
    a = vols(ii,:);
    a = a(~~a);
    
    if length(a)>5
    a = a/sum(a);
    a = sort(a,'descend');
    a = a(1:5);

    axis_vol_score = ( (a - gauss_axis_vols.M) * gauss_axis_vols.C * (a - gauss_axis_vols.M)' );
else
    axis_vol_score = NaN;
    end

    scores(ii) = axis_vol_score;
    
end


