function is_axis_prob = calc_axis_prob(ax_lens)

AX_STD_DEV = 22;
AX_MEAN = 100;

xx = 0:200;

axis_prob = exp(-(xx - AX_MEAN).^2 / (2*AX_STD_DEV^2));

axis_prob = axis_prob / (sqrt(2*pi)*AX_STD_DEV);

ax_lens(ax_lens > max(xx)) = max(xx);

for ii = 1:length(ax_lens)
    
    is_axis_prob(ii) = sum(axis_prob(1:ax_lens(ii)));
    
end