function imshow3(imR, imG, imB, imBG, dim, figNum)

if ~exist('dim', 'var')
    dim = 3;
end

if ~exist('figNum','var')
    figNum = gcf;
end

if length(imR(:))==1
    imR = 0*imB;
elseif length(imG(:))==1
    imG = 0*imB;
elseif length(imB(:))==1
    imB = 0*imG;
end

if exist('imBG') && ~isempty(imBG)
    
    imR = autogain(imR)*.3 + autogain(imBG);
    imG = autogain(imG)*.3 + autogain(imBG);
    imB = autogain(imB)*.3 + autogain(imBG);
    
else
    
    imR = autogain(imR);
    imG = autogain(imG);
    imB = autogain(imB);
    
end

props.dim = dim;
props.sz = size(imR);
props.ind = 1;
props.figNum = figNum;

ff = figure(figNum);clf;
set(ff, 'WindowScrollWheelFcn', @scroller);
set(ff, 'WindowKeyPressFcn', @presser);
% set(ff, 'WindowStyle', 'docked');

drawFrame;


    function drawFrame
        
        if props.dim==1
            imshoww(cat(3, squeeze(imR(props.ind,:,:)),squeeze(imG(props.ind,:,:)),squeeze(imB(props.ind,:,:))));
        elseif props.dim==2
            imshoww(cat(3, squeeze(imR(:,props.ind,:)),squeeze(imG(:,props.ind,:)),squeeze(imB(:,props.ind,:))));
        else
            imshoww(cat(3, squeeze(imR(:,:,props.ind)),squeeze(imG(:,:,props.ind)),squeeze(imB(:,:,props.ind))));
        end
        
        hold on
        text(10,10, num2str(props.ind), 'fontsize', 18, 'color', [1,0,0]);
        drawnow;
        
    end


    function scroller(src,evnt)
        
        props.ind = props.ind + evnt.VerticalScrollCount;
        if props.ind == 0
            props.ind = props.sz(props.dim);
        end

        props.ind = mod(props.ind, props.sz(props.dim));        
        if props.ind == 0
            props.ind = props.sz(props.dim);
        end

        drawFrame;
        
    end


    function presser(src,evnt)
        
        % up and down keys
        props.ind = props.ind + (evnt.Character==30) - (evnt.Character==31);
        
        if props.ind == 0
            props.ind = props.sz(props.dim);
        end
        
        props.ind = mod(props.ind, props.sz(props.dim));
        
        if props.ind == 0
            props.ind = props.sz(props.dim);
        end
            
        if char(evnt.Character) == ' '
            close(props.figNum);
            ff = figure(props.figNum);clf;
            set(ff, 'WindowScrollWheelFcn', @scroller);
            set(ff, 'WindowKeyPressFcn', @presser);
            set(ff, 'WindowStyle', 'docked');
        end
        
        drawFrame;
        
    end

end
