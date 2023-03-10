function nd = analyzeSkel(nd)

% tracing on a single nucleus

crop = nd.traceParams.CROP;

sdata = nd.sdata;

opts = defineSnakeOptsForRAD51OMX;


for jj = 1:length(sdata)
    
    sd = sdata(jj);
    
    sd.snakeTrace = [];
    
    % don't trace short skeletons or those with an intersection
    if sd.length > 10
        
        skk = sd.sk_c;
        mkk = sd.mk;
        
        imm = nd.im(sd.rr, sd.cc, sd.zz);
        
        imm_ = zeros(size(imm)+crop*2);
        skk_ = zeros(size(imm)+crop*2);
        mkk_ = zeros(size(imm)+crop*2);

        imm_(crop+1:end-crop, crop+1:end-crop, crop+1:end-crop) = imm;
        skk_(crop+1:end-crop, crop+1:end-crop, crop+1:end-crop) = skk;
        mkk_(crop+1:end-crop, crop+1:end-crop, crop+1:end-crop) = mkk;

        imm = imm_;
        skk = skk_;
        mkk = mkk_;

        axisData(1).im    = imm;
        axisData(1).skC   = skk;
        axisData(1).mk    = mkk;
        axisData(1).fname = nd.dir;
            
        trace = snakeTrace3v1(axisData, 1, opts);
            

                 
        offset = [min(sd.rr), min(sd.cc), min(sd.zz)] - crop - 1;
        
        if ~isempty(trace) 

            sd.snakeTrace = trace + repmat(offset, size(trace,1), 1);
            
        end
    end
    
    sdd(jj) = sd;
    
end

nd.sdata = sdd;

end





