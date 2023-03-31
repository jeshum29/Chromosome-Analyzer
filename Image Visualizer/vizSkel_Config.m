function vizNucleus(skL, intL, sd, groupSegs, threeWayIntersections, twoWayIntersections, nn, N)


[configArray, unconnectedSegs] = makeConfigArray(threeWayIntersections, nn, N);

% [segsOff, hardFail] = parseConfig(sd, threeWayIntersections, configArray, unconnectedSegs);

segsOff = [];

if ~isempty(twoWayIntersections)
    configArray = cat(1, twoWayIntersections, configArray);
end

groupArray = calcGroupArrayFromConfigArray(configArray, segsOff);

colors = colormap(hsv(128));

colors= ...
    [[0          0.2500    1];...           %blue
    [1           0.1200    0.1200 ];...   %red
    [1           0.7500    0.2000];...     %yellow
    [1           0.2500    1];...            %pink
    [1           0.2000    0];...            %copper
    [0.2000     1           1];...
    [0.2500    0.8000    0.2500];...    %green
    ];           %turquoise


try
    close(11);
end


props.vizVolume = 1;
props.plotTrace = 1;
props.text = 0;

ff = figure(11);clf;
set(ff, 'WindowScrollWheelFcn', @scroller);
set(ff, 'WindowKeyPressFcn', @presser);
drawFrame;


    function drawFrame
        
        figure(1);
        clf;
        
        if props.vizVolume
            
%             mk = bwlabeln(props.nd.sk, 26);
%             mk = ~imdilate(mk==4, ones(3,3,3));
            
            mk = 0*skL + 1;

            % plot the whole skeleton
            vizVolume2(...
                (~~skL) .* mk, ...
                0, ...
                [1,1,1]/2, 1, 0, 1);
            
            
            % intersections and end points mask
            
            mki = ~((~~intL) + (count_neighbor_pixels(double(~~skL), 26) == 1));
            mki = mk;
            
            hold on
            
            sz = size(groupArray);
            segsSingle = [];
            
            % count single segments
            for ii = 1:length(groupSegs)
                groupSegs(ii)
                if isempty(intersect(groupSegs(ii), groupArray(:))) &&...
                        sd(groupSegs(ii)).flag &&...
                        isempty(intersect(groupSegs(ii), segsOff))
                    
                    segsSingle = [segsSingle, groupSegs(ii)];
                end
            end
            
            groupSegs = setxor(segsSingle, groupSegs);
            groupArrayTmp = groupArray;
            counter = 0;
            
            % count connected groups
            for ii = 1:length(groupSegs)
                [row, col] = ind2sub(sz, find(groupArrayTmp==groupSegs(ii)));
                if isempty(row), continue, end
                counter = counter+1;
                groupArrayTmp(row, :) = 0;
            end
            
            % get a list of colors
            numberOfColors = counter + length(segsSingle);
            interval = floor(size(colors,1)/(numberOfColors+1));
            cc = colors(1:interval:end,:);
            
            % plot isolated segments in the group
            for ii = 1:length(segsSingle)
                
                vizVolume2(...
                    (skL==segsSingle(ii) ).*mki,...
                    0,...
                    cc(ii,:), 1, 0, 1/2);
                
            end
            
            groupArrayTmp = groupArray;
            counter = ii;
            
            if isempty(counter)
                counter = 0;
            end
            
            % plot connected segments in the same color
            for ii = 1:length(groupSegs)
                [row, col] = ind2sub(sz, find(groupArrayTmp==groupSegs(ii)));
                if isempty(row), continue, end
                counter = counter+1;
                subsegs = unique(nonzeros(groupArrayTmp(row, :)));
                
                for jj = 1:length(subsegs)
                    if isempty(intersect(subsegs(jj), segsOff))
                        
                        vizVolume2(...
                            (skL==subsegs(jj)).*mki,...
                            0,...
                            cc(counter,:), 1, 0, 1/2);
                        
                    end
                end
                
                groupArrayTmp(row, :) = 0;
            end
            
            
            % intersections and end points in white
            vizVolume2(...
                ~~intL,...
                0,...
                [1,1,1], 1, 0, 1);
            
            
        end
        
        
        if props.text
            
            for ii = 1:length(sd)
                
                text(mean(sd(ii).cc)+2, mean(sd(ii).rr)+2, mean(sd(ii).zz)+2, num2str(round(max(sd(ii).intensities))),'color', 'w', 'fontsize', 16);
                
            end
        end
        
        
        if props.plotTrace
            
            for jj = 1:length(sd)
                
                if ~isfield(sd(jj), 'trace')
                    continue;
                end
                
                if isempty(sd(jj).trace)
                    continue;
                end
                
                plot3(...
                    sd(jj).trace(:,2), ...
                    sd(jj).trace(:,1), ...
                    sd(jj).trace(:,3), ...
                    '-', 'LineWidth', 5, 'color', colors(1,:));
                
                axis equal
                axis vis3d
                
            end
        end
        %   text(10,10, num2str(props.ind), 'fontsize', 18, 'color', [1,0,0]);
        
        drawnow;
        
    end


    function scroller(src,evnt)
        
        if evnt.VerticalScrollCount > 0
            props.ind = min(props.ind + 1, props.sz);
        elseif evnt.VerticalScrollCount < 0
            props.ind = max(props.ind - 1, 1);
        end
        drawFrame;
        
    end


    function presser(src,event)
        
        if char(event.Character) == 'v'
            props.vizVolume = ~props.vizVolume;
        elseif char(event.Character) == 'p'
            props.plotTrace = ~props.plotTrace;
        elseif char(event.Character)=='t'
            props.text = ~props.text;
        end
        
        if event.Character == 32
            close(1);
            ff = figure(1);clf;
            set(ff, 'WindowScrollWheelFcn', @scroller);
            set(ff, 'WindowKeyPressFcn', @presser);
        end
        
        drawFrame;
        
    end

end
