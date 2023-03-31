function nd = analyzeSkel(nd)

% tracing on a single nucleus

crop = nd.traceParams.CROP;

% nd = makeSDataSIM(nd);
sdata = nd.sdata;


for jj = 1:length(sdata)
    
    sd = sdata(jj);
    sd.trace = [];
    sd.trace_ = [];
    sd.traceInit = [];
    sd.trace_fail = 0;
    
    % don't trace short skeletons or those with an intersection
    if sd.length > 4
        
        skk = sd.sk_c;
        imm = nd.im(sd.rr, sd.cc, sd.zz);
                
         try
            
        [xff, yff, zff, xf, yf, zf, xi, yi, zi] = ...
                        traceChr3v2(imm, skk, ...
                        nd.traceParams.CROP,...
                        nd.traceParams.INTERPS_PER_PIXEL,...
                        nd.traceParams.STEPS_PER_PIXEL,...
                        nd.traceParams.CONV_LENGTH);
            
         catch
             sd.trace_fail = 1;
             'failed trace!'
         end
                 
        offset = [min(sd.rr), min(sd.cc), min(sd.zz)] - crop - 1;
        
        if ~sd.trace_fail && ~isempty(xf) 
            
            mk = xff~=0;
            sd.trace =  [xff(mk)' + offset(1), yff(mk)' + offset(2), zff(mk)' + offset(3)];
            
            mk = xf~=0;
            sd.trace_ = [xf(mk)' + offset(1), yf(mk)' + offset(2), zf(mk)' + offset(3)];
            
            sd.traceInit = [xi' + offset(1), yi' + offset(2), zi' + offset(3)];
            
        end
    end
    
    sdd(jj) = sd;
    
end

nd.sdata = sdd;

end





