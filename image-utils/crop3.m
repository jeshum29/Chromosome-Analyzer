function imshow3(im, dim)

global imc

if ~exist('dim', 'var')
    dim = 3;
end

w = whos('im');

if strcmp(w.class, 'logical')
    im = double(im);
end

props.im = im;
clear im
props.dim = dim;
props.sz = size(props.im);
props.ind = 1;

ff = figure(1);clf;
set(ff, 'WindowScrollWheelFcn', @scroller);
set(ff, 'WindowKeyPressFcn', @presser);
drawFrame;


    function drawFrame
        
        maxx = max(props.im(:));
        minn = min(props.im(:));
        
        if props.dim==1
            imshow(squeeze(props.im(props.ind,:,:)),[minn,maxx]);
        elseif props.dim==2
            imshow(squeeze(props.im(:,props.ind,:)),[minn,maxx]);
        else
            imshow(squeeze(props.im(:,:,props.ind)),[minn,maxx]);
        end
        
        hold on
        text(10,10, num2str(props.ind), 'fontsize', 18, 'color', [1,0,0]);
                
        
        
        drawnow;
        
    end


    function scroller(src,evnt)
        
        if evnt.VerticalScrollCount > 0
            props.ind = min(props.ind + 1, props.sz(props.dim));
        elseif evnt.VerticalScrollCount < 0
            props.ind = max(props.ind - 1, 1);
        end
        drawFrame;
        
    end


    function presser(src,evnt)
        
        if strcmp(evnt.Key, 'uparrow')
            props.ind = min(props.ind + 1, props.sz(props.dim));
        elseif strcmp(evnt.Key, 'downarrow')
            props.ind = max(props.ind - 1, 1);
        elseif strcmp(evnt.Key, 'space')
            close(1);
            ff = figure(1);clf;
            set(ff, 'WindowScrollWheelFcn', @scroller);
            set(ff, 'WindowKeyPressFcn', @presser);
        elseif char(evnt.Character)=='c'
            [props.y, props.x] = ginput(1);
            props.z = props.ind;
        elseif char(evnt.Character)=='v'
            [props.y2, props.x2] = ginput(1);
            props.z2 = props.ind;
        elseif strcmp(evnt.Key, 'return')
            sz = props.sz;
            y = min(props.y, props.y2);
            x = min(props.x, props.x2);
            y2 = max(props.y, props.y2);
            x2 = max(props.x, props.x2);
            z = min(props.z, props.z2);
            z2 = max(props.z, props.z2);
            imc = props.im(max(1,round(x)):min(sz(1),round(x2)), max(1,round(y)):min(sz(2),round(y2)), z:z2);
        end
        
        drawFrame;
        
    end

end
