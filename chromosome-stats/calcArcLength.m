function arcLength = calcArcLength(trace)


steps = trace(2:end,:) - trace(1:end-1,:);
steps = sqrt(sum(steps.^2,2));

arcLength(1) = 0;

for ii = 1:size(steps,1)
    
    arcLength(ii+1) = sum(steps(1:ii));
    
end