function sk = deleteSpurs(sk, num)

for ii = 1:num
    
    spurs = (count_neighbor_pixels(double(sk), 26) == 1);
    
    sk(spurs) = 0;
    sk = skel3(double(sk));

end

