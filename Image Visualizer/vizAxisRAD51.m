function vizAxesRAD51(axisData, ind, max_spotAxisDistance, min_spotScore)

figNum = 2718;

figure(figNum);
clf;

ii = ind;

axisData(ii).fname

fname = axisData(ii).fname;
fname = fname(1: (regexp(fname, 'nukes') - 1));

gonad = load([fname 'gonad.mat']);
gonad = gonad.gonad;


imRGB = makeAxisImageFISH(axisData, ind);

vizIm(imRGB, figNum, []);
set(gcf, 'Alphamap', (0:64)*.3/64);


hold on

vizVolume3(axisData(ii).skC, [250  120  30]/255, figNum, 0, 1, .3); %, min(cc), min(rr), min(zz))


try
    traceColor = [0,1,0];
    plot3(...
        axisData(ii).snakeTrace(:,2),...
        axisData(ii).snakeTrace(:,1),...
        axisData(ii).snakeTrace(:,3),...
        '-o',...
        'LineWidth', 2,...
        'MarkerSize', 3,...
        'MarkerFaceColor', traceColor,...
        'MarkerEdgeColor', traceColor,...
        'color', traceColor);
end


spots = [];
if isfield(axisData, 'spots_im1')
    if ~isempty(axisData(ii).spots_im1)
        spots = axisData(ii).spots_im1;
    end
end

if isfield(axisData, 'spots_im2')
    if ~isempty(axisData(ii).spots_im2)
        spots = axisData(ii).spots_im2;
    end
end

spot_offset = -1/2;
spot_scale = 1;

rr = axisData(ii).rr;
cc = axisData(ii).cc;
zz = axisData(ii).zz;

for jj = 1:length(spots)
    
    goodSpotFlag = 0;
    plotSpotFlag = 0;
    
    if ~isempty(axisData(ii).spotAxisDistances)  
        plotSpotFlag = axisData(ii).spotAxisDistances(jj) < 15 && spots(jj).intensity_score > 10;
        goodSpotFlag = axisData(ii).spotAxisDistances(jj) < max_spotAxisDistance && spots(jj).intensity_score > min_spotScore;
    end
    
    if ~plotSpotFlag
        continue;
    end
    
    if goodSpotFlag
        spotColor = [0,1,0];
    else
        spotColor = [1,1,0];
    end
    
    if plotSpotFlag
        
        pos = [...
                spots(jj).r(2) * spot_scale - min(cc) + spot_offset,...
                spots(jj).r(1) * spot_scale - min(rr) + spot_offset,...
                spots(jj).r(3) * spot_scale - min(zz) + spot_offset,...
              ];
        
        scatter3(pos(1), pos(2), pos(3), 200, spotColor, 'o', 'filled');
        
        text_str = [...
                    num2str(round(spots(jj).intensity_score*10)/10), ', ',...
                    num2str(round(axisData(ii).spotAxisDistances(jj)*10)/10)];
        
        text_offset = 3;
        
        text(pos(1) + text_offset, pos(2) + text_offset, pos(3) + text_offset, text_str, 'fontsize', 12, 'color', spotColor);
        
    end
    
    
    '';
    
    
end


