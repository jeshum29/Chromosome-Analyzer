function axisData = batchCatAxes(dds)

% code for extracting good axes from OMX images of RAD51 IF
% Mostly the same as the original code for DV images

% Sept 2013

slash = filesep;

axisData = [];

   for nn = 1:length(dds)
        
       if isempty(dds(nn).name)
           continue;
       end
       
        gonad = load([dds(nn).name slash 'gonad.mat']);
        gonad = gonad.gonad;
        
        if ~isfield(gonad, 'finalized_nukes')
            continue;
        end
        
        dds(nn).name
        
        for ii = find(gonad.finalized_nukes==1)
            
           fname = [ gonad.writeDir filesep ...
                     'nukes' filesep ...
                     num2str(gonad.good_nukes(ii)) filesep ...
                     'nd.mat']
                 
           fname = [ dds(nn).name filesep ...
                     'nukes' filesep ...
                     num2str(gonad.good_nukes(ii)) filesep ...
                     'nd.mat']
                 
           % these are some of the nuclei that cleanSkel2 gets stuck on
           % (presumably there are loops present)
           if strcmp(fname, 'E:\Data\2014_11_12_syIs44_RAD51_realign\01\2014_11_12_syIs44_RAD51_01_02_SIR_ALX\nukes\42\nd.mat')
               continue;
           end
           
           if strcmp(fname, 'E:\Data\2014_11_12_syIs44_RAD51_realign\02\2014_11_12_syIs44_RAD51_02_01_SIR_ALX\nukes\22\nd.mat')
               continue;
           end
            
            nd = load(fname);
            nd = nd.nd;
            
            if ~isfield(nd, 'mkL_final_includeFlags')
                continue;
            end
            
            try
                ad = makeAxisDataForFISH(nd);
            catch
                'ERROR in makeAxisDataForFISH'
                ad = [];
            end
            
            axisData = [axisData, ad];
                   
            
                   
        end
   end
   
   '';
   
        
        