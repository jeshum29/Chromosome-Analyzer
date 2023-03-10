function batchWriteTIFFfromDV(rootDir, writeProjection, writeStack, FILTER_STRING)


% if DVs are in top-level/root directory:
dvFiles = dir([rootDir filesep '*' FILTER_STRING '*.dv']);

for ii = 1:length(dvFiles)
    dvFiles(ii).name = [rootDir filesep dvFiles(ii).name];
end

% load each dv file and write either a TIFF stack or a projection 
for ii = 1:length(dvFiles)

    if dvFiles(ii).isdir
        continue;
    end
        
    writeTIFFfromOMXDV([dvFiles(ii).name], writeProjection, writeStack);

end

