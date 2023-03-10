function nd = vizAxisTrace(nd, figNum, clfFlag, cc)

if ~exist('figNum', 'var')
    figNum = 99;
end

if ~exist('clfFlag','var')
    clfFlag = 1;
end

if ~exist('cc','var')
    cc = ones(length(nd.sdata), 3);
end

if size(cc, 1)==1
    cc = repmat(cc, length(nd.sdata), 1);
end

ff = figure(figNum);

if clfFlag
    clf;
end

set(ff, 'color', [0,0,0]);


for jj = 1:length(nd.sdata)
    
    if ~isfield(nd.sdata(jj), 'trace') || isempty(nd.sdata(jj).trace)
        continue;
    end
    
    hold on
    
    plot3(...
        nd.sdata(jj).trace(:,2), ...
        nd.sdata(jj).trace(:,1), ...
        nd.sdata(jj).trace(:,3), ...
        '-', 'LineWidth', 3,...
        'MarkerSize', 1,...
        'MarkerFaceColor', cc(nd.sdata(jj).label,:),...
        'MarkerEdgeColor', cc(nd.sdata(jj).label,:),...
        'color', cc(nd.sdata(jj).label,:));
    
%     plot3(...
%         nd.sdata(jj).traceInit(:,2), ...
%         nd.sdata(jj).traceInit(:,1), ...
%         nd.sdata(jj).traceInit(:,3), ...
%         '-o', 'LineWidth', 2,...
%         'MarkerSize', 1,...
%         'MarkerFaceColor', cc(nd.sdata(jj).label,:),...
%         'MarkerEdgeColor', cc(nd.sdata(jj).label,:),...
%         'color', cc(nd.sdata(jj).label,:)/3);
end


daspect([1 1 1]);
view(3);
axis tight;
camlight left;
camlight right;
axis vis3d;
lighting none; %without this, no intensity mapping
%lighting flat;
material dull;

rotate3d;

set(gca, 'Position', [0,0,1,1]);

% 
% set(gca, 'color', 'w');
% set(gcf, 'color', 'w');

axis off



