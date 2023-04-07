function batchMakeNukesFromMask4(dds)

for ii = 1:length(dds)
    
    dds(ii).name
    
    gonad = load([dds(ii).name filesep 'gonad.mat']);
    gonad = gonad.gonad;
    
    bandpass = load3([dds(ii).name filesep 'mk_bandpass.tif']);
    bandpass = bandpass;

    makeNukesFromMask4(bandpass, gonad);
    
end