function batchInitializeGonadOMX(dds, NUM_CHAN, SEG_CHAN, OVERWRITE_FLAG)

if ~exist('OVERWRITE_FLAG', 'var')
    OVERWRITE_FLAG = 0;
end



for ii = 1:length(dds)
    
    initializeGonadOMX(dds(ii).name, NUM_CHAN, SEG_CHAN, OVERWRITE_FLAG);
   
    
end