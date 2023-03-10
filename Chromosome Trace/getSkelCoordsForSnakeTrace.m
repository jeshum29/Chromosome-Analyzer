

function [pixelCoords] = getSkelCoords_internal(sk)

x = []; y = []; z = [];

pixelCoords = [];

numNeighbors = find_intersections_3(double(sk), 26);
start = find(numNeighbors(:)==1,1,'first');

if isempty(start)
    return;
end

[x_s, y_s, z_s] = ind2sub(size(sk), start);

x(1) = x_s;
y(1) = y_s;
z(1) = z_s;

sk(x(end),y(end),z(end)) = 0;

try
    nhood = sk(x(end)-1:x(end)+1,y(end)-1:y(end)+1,z(end)-1:z(end)+1);
catch
    x = [];
    y = [];
    z = [];
    return;
end

go = 1;

while go
    
    next = find(nhood,1,'first');
    [x_s, y_s, z_s] = ind2sub([3,3,3], next);
    
    x(end+1) = x(end) + x_s - 2;
    y(end+1) = y(end) + y_s - 2;
    z(end+1) = z(end) + z_s - 2;
    
    sk(x(end),y(end),z(end)) = 0;
    
    % if this line fails, the skeleton borders the image
    try
        nhood = sk(x(end)-1:x(end)+1,y(end)-1:y(end)+1,z(end)-1:z(end)+1);
    catch
        x = [];
        y = [];
        z = [];
        break;
    end
    
    if numNeighbors(x(end),y(end),z(end))==1
        go = 0;
    end
    
end

pixelCoords = [x', y', z'];

end

