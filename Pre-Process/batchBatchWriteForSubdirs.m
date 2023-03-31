function dvFileDirs = batchBatchWriteForSubdirs(rootDir, FILTER_STRING, WRITE_FILES_FLAG)

% function to run batchWriteTIFFfromOMXDV, etc, on the DV files
% in each subdirectory (usually corresponding to a single gonad)
% within rootDir

subDirs = dir([rootDir filesep]);
dvFileDirs = {}; 
count = 0;

for ii = 3:length(subDirs)
    if isdir([rootDir filesep subDirs(ii).name])

        subDirFullPath = [rootDir filesep subDirs(ii).name filesep];
        
        if WRITE_FILES_FLAG
            
            tmp = dir([subDirFullPath filesep '*.tif']);
            if length(tmp)
                % continue;
            end
            
            % write projections and stacks
            batchWriteTIFFfromOMX(subDirFullPath, 1, 1, FILTER_STRING);

            % write isotropic TIFF stacks
            batchWriteIsotropicTIFFfromOMXTIFF(subDirFullPath, [], FILTER_STRING);

        end
        
        % here we generate a list of DV file subdirectories -- used for
        % initializeGonadOMX 
        
        dvFileDirs_ = dir([rootDir filesep subDirs(ii).name filesep])';

        % append the complete path
        for jj = 3:length(dvFileDirs_)
            if dvFileDirs_(jj).isdir
                
                count = count + 1;
                dvFileDirs(count).name = [rootDir filesep subDirs(ii).name filesep dvFileDirs_(jj).name];
                
                disp(['dds(' num2str(count) ').name = ''' dvFileDirs(count).name ''';'])
                
            end
        end
    end
end



