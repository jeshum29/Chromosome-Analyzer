function batchMakeBackgroundMask(dds)

for ii = 1:length(dds)
    
    dds(ii).name
    
    gonad = load([dds(ii).name filesep 'gonad.mat']);
    gonad = gonad.gonad;
    
    makeBackgroundMask(gonad);
    
end

