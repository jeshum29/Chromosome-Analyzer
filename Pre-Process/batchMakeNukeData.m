function batch(dds, OVERWRITE_FLAG)

N_total = 0;

for ii = 1:length(dds)
    
    ii
    
    if isempty(dds(ii).name)
        continue;
    end
    
    try
        gonad = load([dds(ii).name filesep 'gonad.mat']);
    catch
        continue;
    end
    
    gonad = gonad.gonad;
    
%     gonad.NUM_CHAN = 2;
%     save([gonad.writeDir filesep 'gonad.mat'], 'gonad');

    crop = 5;
    
    %try
    N = makeNukeData(gonad, crop, OVERWRITE_FLAG);
    %catch
     %   'WARNING: makeNukeData error'
      %  dds(ii).name
       % N = 0;
    end
    
    N_total = N + N_total
    
end

