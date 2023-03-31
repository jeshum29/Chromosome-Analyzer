function mov = imshow3(im, dim, speed)

if ~exist('dim', 'var')
    dim = 3;
end

if ~exist('speed','var')
    speed = 1;
end

w = whos('im');

if strcmp(w.class, 'logical')
    im = double(im);
end

figure;


    maxx = max(im(:));
    minn = min(im(:));

    for ii = 1:size(im,dim)

        if dim==1
            imshow(squeeze(im(ii,:,:)),[minn,maxx]);
        elseif dim==2
            imshow(squeeze(im(:,ii,:)),[minn,maxx]);
        else
            imshow(squeeze(im(:,:,ii)),[minn,maxx]);
        end

        
        hold on
        
        text(10,10, num2str(ii), 'fontsize', 18, 'color', [1,0,0]);
        
        drawnow;
        mov(ii) = getframe;
        for jj = 1:10^speed; end

       '';
    end



