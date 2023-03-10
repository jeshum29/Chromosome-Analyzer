
function subGroup = makeSubGroups(regRegArray, regs)

% returns a structure with an entry for each set of good connected segments
% (bad segments get an NaN in regRegArray, which this function never sees)

ssa = regRegArray;
try
isolatedRegs = intersect(find(ssa(:,1)==0), regs);
catch
    ''
end

if ~isempty(isolatedRegs)
    
    for ii = 1:length(isolatedRegs)
        subGroup(ii).regs = isolatedRegs(ii);
    end
    subGroupCounter = length(isolatedRegs);
    
else
    subGroupCounter = 0;
end


go_nextSubGroup = 1;
while go_nextSubGroup
    
    subGroupCounter = subGroupCounter + 1;
    subGroup(subGroupCounter).regs = [];
    
    current_regs = find(ssa(:,1)>0,1,'first');
    
    go_traceSubGroup = 1;
    while go_traceSubGroup
        
        neighbors = unique(nonzeros(ssa(current_regs, :)));
        ssa(current_regs, :) = 0;
        
        subGroup(subGroupCounter).regs = unique([subGroup(subGroupCounter).regs, neighbors']);
        current_regs = neighbors;
        
        if isempty(current_regs)
            go_traceSubGroup = 0;
        end
        
    end
    
    
    if ~sum(ssa(~isnan(ssa)))
        go_nextSubGroup = 0;
    end
    
end

