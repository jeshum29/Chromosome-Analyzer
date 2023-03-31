function spot = spotter3(im, mk_axis, thresh, smoothing)

% ----------------------------------------------------
%
% fit 3D gaussian to the foci in the image
%
% ----------------------------------------------------


slash = filesep;

opt =  optimset('MaxIter', 10, 'Display', 'off',  'TolX',  1/10, 'TolFun',  1e-8);

INIT_SPOT_WIDTH = 2*  [1, 1, 1];
MIN_SPOT_WIDTH  = 1/4*[1, 1, 1];
MAX_SPOT_WIDTH  = 6*  [1, 1, 1];

CROP   = [4, 4, 8];
SBOUND = [2, 2, 4];

disp_flag = 1;
spot_num  = 0;

im_orig = im;

im_pad = zeros(size(im) + CROP*2);
mk_pad = im_pad;

im_pad( CROP(1)+1 : end-CROP(1), CROP(2)+1 : end-CROP(2), CROP(3)+1 : end-CROP(3) ) = im;
mk_pad( CROP(1)+1 : end-CROP(1), CROP(2)+1 : end-CROP(2), CROP(3)+1 : end-CROP(3) ) = mk_axis;

im =      im_pad;
mk_axis = mk_pad;

im = double(im);

im_back  = im - mean(im(:));
im       = double((im));
I_bar    = mean(im(:));
I_std    = std(im(:));

mk_back = imerode(im_back~=0, ones(2, 2, 2));

mk_back = mk_back .* mk_axis;

imF = imfilter3(im, [1,1,1]*smoothing);

ws = double(watershed( -imF.*mk_back, 26) ) .*(mk_back);

pp = regionprops(ws, 'Area');

for ii = 1:length(pp)
    
    maxx   = max(im(ws(:)==ii));
    zScore = (maxx - I_bar)/I_std;
    
    if zScore < thresh
        ws(ws==ii) = 0;
    else        
        sum(im(ws(:)==ii) > maxx/2);
    end
end

mk  = ~~ws;
mkL = bwlabeln(mk);
sz  = size(mkL);
num_segs = max(mkL(:));

if 0
    
    imshow3ck(autogain(~mkL)/3, 0, 0, autogain(im), 3, 33);
    
    '';
end

for jj = 1:num_segs
    
    mkk     = (mkL==jj).*(im~=0);
    imF_mk  = imF.*mkk;
    
    [maxValue, ind] = max(imF_mk(:));
    [x, y, z]       = ind2sub(sz, ind);
    
    INIT_INTENSITY  = im(x, y, z);
    MIN_INTENSITY   = 0;
    MAX_INTENSITY   = INIT_INTENSITY*2;
    
    INIT_OFFSET     = I_bar;
    MIN_OFFSET      = max( 0, INIT_OFFSET - 3*std(im(~~mkk(:))) );
    MAX_OFFSET      = INIT_OFFSET + 3*std(im(~~mkk(:)));
    
    INIT_CENTER     = [x, y, z];
    
    ddd  = [ INIT_CENTER,          INIT_INTENSITY, INIT_SPOT_WIDTH, INIT_OFFSET ];
    lddd = [ INIT_CENTER - SBOUND, 0,              MIN_SPOT_WIDTH,  MIN_OFFSET  ];
    uddd = [ INIT_CENTER + SBOUND, MAX_INTENSITY,  MAX_SPOT_WIDTH,  MAX_OFFSET  ];
    
    rr = x + ( -CROP(1):CROP(1) );
    cc = y + ( -CROP(2):CROP(2) );
    zz = z + ( -CROP(3):CROP(3) );
    
    [rrm, ccm, zzm] = ndgrid(rr, cc, zz);
    
    mkkk = mkk(rr, cc, zz);
    imm  = im(rr, cc, zz);
    
    %     mkkk = mkkk.*(imm > INIT_INTENSITY/4);
    
    fit_init = g_fit_fun(rrm, ccm, zzm, ddd);
    
    errorFlag = 0;
    try
        [ddd, resn, res, exit_flag, output, lambda, jac] = lsqnonlin( @fitter, ddd, lddd, uddd, opt );
    catch
        errorFlag = 1;
        'Warning: lsqnonlin error!'
    end
    
    fit_final = g_fit_fun(rrm, ccm, zzm, ddd);
    
    xp = find(rr==x);
    yp = find(cc==y);
    zp = find(zz==z);
    
    if 0
        
        %blue, red, green
        %
        %         figure(10);clf;
        %         plot([squeeze(imm(yp,xp,:)), squeeze(fit_init(yp,xp,:)), squeeze(fit_final(yp,xp,:))]);
        %         title('Z');
        %         figure(11);clf;
        %         plot([squeeze(imm(yp,:,zp))', squeeze(fit_init(yp,:,zp))', squeeze(fit_final(yp,:,zp))']);
        %         title('Y');
        %         figure(12);clf;
        %         plot([squeeze(imm(:,xp,zp)), squeeze(fit_init(:,xp,zp)), squeeze(fit_final(:,xp,zp))]);
        %         title('X');
        
        imm((imm - I_bar - 3*I_std) < 0) = 0;
        vizIm(imm,99);
        colormap(gray(255));
        set(gcf, 'Alphamap', (0:64)*(.5)/64);
        hold on
        
        offset = 1/2;
        r = ddd(1:3) - [rr(1) cc(1) zz(1)] + offset;
        scatter3( r(2), r(1), r(3), 100, [1,0,0], 'o', 'filled');
        
        %         im_viz = im - I_bar - 3*I_std;
        %         im_viz(im_viz < 0) = 0;
        %         vizIm(autogain(im_viz), 99);
        %         colormap(gray(255));
        %         hold on
        %         scatter3(ddd(1) - 1/2, ddd(2) - 1/2, ddd(3) - 1/2, 100, [1,0,0], 'o', 'filled');
        
        %     proj = max(im_viz,[],3);
        %     figure(99);clf;
        %     imshow(autogain(proj),[],'InitialMagnification', 1000);
        %     hold on
        %     scatter(ddd(1), ddd(2), 30, [1,0,0], 'o', 'filled');
        
        
        %         im_disp = imm;
        %
        %         figure(1000);clf;
        %         projj = squeeze(max(im_disp.*mkkk,[],3));
        %         imagesc(cc,rr,autogain(projj));
        %         colormap(gray);
        %         hold on
        %         scatter(ddd(2),ddd(1),30,[1,0,0],'o','filled');
        %
        %         figure(1110);clf;
        %         projj = squeeze(max(im_disp.*mkkk,[],1));
        %         imagesc(zz,cc,autogain(projj));
        %         colormap(gray);
        %         hold on
        %         scatter(ddd(3),ddd(2),30,[1,0,0],'o','filled');
        %
        %         figure(1220);clf;
        %         projj = squeeze(max(im_disp.*mkkk,[],2));
        %         imagesc(zz,rr,autogain(projj));
        %         colormap(gray);
        %         hold on
        %         scatter(ddd(3),ddd(1),30,[1,0,0],'o','filled');
        
        '';
    end
    
    xpos = ddd(1);
    ypos = ddd(2);
    zpos = ddd(3);
    
    final_intensity   = ddd(4);
    final_spot_width  = ddd(5:7);
    final_offset      = ddd(8);
    
    res_score       = resn;
    intensity_score = (final_intensity - I_bar)/I_std;
    
    fit_final = fit_final.*mkkk;
    
    Itotal = sum(double((fit_final(:)-final_offset).*mkkk(:)));
    % Ivar = I_std*sqrt(sum(double(logical(mkkk(:)))));
    
    Iwidth         = sqrt(sum(ddd(5:7).^2));
    integral_score = Itotal/Iwidth;
    
    '';
    
    if intensity_score > 0 && ~errorFlag
        
        spot_num                        = spot_num+1;
        spot(spot_num).r                = ddd(1:3) - CROP;
        spot(spot_num).intensity        = final_intensity;
        spot(spot_num).offset           = final_offset;
        spot(spot_num).intensity_score  = intensity_score;
        spot(spot_num).integral_score   = integral_score;
        spot(spot_num).res_score        = res_score;
        spot(spot_num).width            = final_spot_width;
        
    end
end

if spot_num==0
    spot = [];
end

if 0
    
    mp = max(im,[],3);
    figure(3);clf;
    imshow(autogain(mp),[],'initialmagnification', 1000);
    
    hold on
    
    for jj = 1:length(spot)
        scatter(spot(jj).r(1), spot(jj).r(2),30,[1,0,0],'o','filled');
    end
end

'';

if 0

    im_viz = im_orig - I_bar - 3*I_std;
    im_viz(im_viz < 0) = 0;
    
    vizIm(autogain(im_viz), 88);
    cc = colormap(gray(255));
    colormap(cc);
    set(gcf, 'Alphamap', (0:64)*(.5)/64);
    
    hold on
    offset = -1/2;
    
    nuke_cen = mean(reshape([spot.r],3,length([spot.r])/3)', 1);
    
    cc = [.2,1,.2];
    
    for jj = 1:length(spot)
        
        
        offset = -1/2*[1,1,1];
        scatter3(spot(jj).r(2)+offset(2), spot(jj).r(1)+offset(1), spot(jj).r(3)+offset(3), 100, cc, 'o', 'filled');
        
        offset = 12;
        unit_vec = spot(jj).r - nuke_cen;
        unit_vec = unit_vec ./ sqrt(sum(unit_vec.*unit_vec));
        
        text_pos = spot(jj).r + offset * unit_vec;
        
        text_str = [...
            num2str(round(spot(jj).intensity_score*100)/100), ', ',...
            num2str(round(spot(jj).intensity/1e3)), ', ',...
            num2str(round(spot(jj).integral_score*100)/100), ', ',...
            ];
        
        text_str = [...
            num2str(round(spot(jj).intensity_score*100)/100), ', ',...
            ];
        
        text(...
            text_pos(2), ...
            text_pos(1), ...
            text_pos(3), ...
            text_str, 'fontsize', 12, 'color', cc);
        
        plot3(...
            [text_pos(2),  spot(jj).r(2)],...
            [text_pos(1),  spot(jj).r(1)],...
            [text_pos(3),  spot(jj).r(3)],...
            'color', cc);

    end
    
    '';
end


    function C = fitter(ddd)
        
        gg = g_fit_fun(rrm,ccm,zzm,ddd);
        
        tmp = (double(imm)-gg);
        C = tmp(logical(mkkk));
        
    end

end


function C = g_fit_fun(X,Y,Z,ddd)

x0 = ddd(1); 
y0 = ddd(2); 
z0 = ddd(3);

bx = ddd(5); 
by = ddd(6); 
bz = ddd(7);

IG = ddd(4);
I0 = ddd(8);

C = I0 + IG*exp( -(...
                    (X - x0).^2 / (2*bx^2) + ...
                    (Y - y0).^2 / (2*by^2) + ...
                    (Z - z0).^2 / (2*bz^2)...
                    ) );

end

% function C = g_fit_fun(X,Y,x0,y0,IG,I0,b1,theta,b2)
%
% rot = [[cos(theta), sin(theta)];[-sin(theta), cos(theta)]];
% KK = rot*diag([b1^-2,b2^-2])*rot'/2;
% K11 = KK(1,1);
% K12 = KK(1,2);
% K22 = KK(2,2);
%
% C = I0 + IG*exp( -(K11*(X-x0).^2+2*K12*(Y-y0).*(X-x0)+K22*(Y-y0).^2) );
%
% end

