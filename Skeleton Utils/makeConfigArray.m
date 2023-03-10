function [configArray, unconnectedSegs] = makeConfigArray(intersections, nn, N)

        config = sprintf(['%0' num2str(N) 'd'], str2num(dec2base(nn, 4)));
        configArray = [];
        
        for jj = 1:N
            
            cc = str2num(config(jj));
            switch cc
                
                case 0
                    configArray(jj,1:2) = 0;
                    unconnectedSegs(jj) = 0;
                case 1
                    configArray(jj,1:2) = intersections(jj, [1,2]);
                    unconnectedSegs(jj) = intersections(jj, 3);
                case 2
                    configArray(jj,1:2) = intersections(jj, [2,3]);
                    unconnectedSegs(jj) = intersections(jj, 1);                    
                case 3
                    configArray(jj,1:2) = intersections(jj, [1,3]);
                    unconnectedSegs(jj) = intersections(jj, 2);                    
                    
            end
        end