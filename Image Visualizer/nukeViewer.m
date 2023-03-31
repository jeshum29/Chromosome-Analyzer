function varargout = nukeViewer(varargin)
% NUKEVIEWER MATLAB code for nukeViewer.fig
%      NUKEVIEWER, by itself, creates a new NUKEVIEWER or raises the existing
%      singleton*.
%
%      H = NUKEVIEWER returns the handle to a new NUKEVIEWER or the handle to
%      the existing singleton*.
%
%      NUKEVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NUKEVIEWER.M with the given input arguments.
%
%      NUKEVIEWER('Property','Value',...) creates a new NUKEVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before nukeViewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to nukeViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help nukeViewer
% Last Modified by GUIDE v2.5 04-Jun-2014 15:35:50
% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @nukeViewer_OpeningFcn, ...
    'gui_OutputFcn',  @nukeViewer_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);

if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

end

% --- Executes just before nukeViewer is made visible.
function nukeViewer_OpeningFcn(hObject, eventdata, handles, varargin)

if length(varargin)==1
    handles.dd = varargin{1};
else
    handles.dd = get(handles.textBox_dir,'String');
end

handles.gonad = load([handles.dd filesep 'gonad.mat']);
handles.gonad = handles.gonad.gonad;

set(handles.textBox_dir, 'String', handles.dd);

handles = load_gonad(handles);
set(handles.listbox_surfs, 'String', handles.nukeFoldersList, 'Value', 1);

display_gonad(handles, hObject);

% Update handles structure
handles.output = hObject;
guidata(hObject, handles);


end

% --- Outputs from this function are returned to the command line.
function varargout = nukeViewer_OutputFcn(hObject, eventdata, handles)
% Get default command line output from handles structure
varargout{1} = handles.output;

end


function textBox_dir_Callback(hObject, eventdata, handles)

handles.dd     = get(hObject,'String');
handles.gonad  = load([handles.dd filesep 'gonad.mat']);
handles.gonad  = handles.gonad.gonad;

handles = load_gonad(handles);

set(handles.listbox_surfs, 'String', handles.nukeFoldersList, 'Value', 1);

display_gonad(handles, hObject);
handles.index_selected = 1;
guidata(hObject, handles);

end

% =========================================================-
%
%  INITIALIZATION
%
% =========================================================-


function handles = load_gonad(handles)

try
    set(handles.im1_scale,'String',num2str(handles.gonad.im1_scale));
    set(handles.im2_scale,'String',num2str(handles.gonad.im2_scale));
    set(handles.im3_scale,'String',num2str(handles.gonad.im3_scale));
    set(handles.im4_scale,'String',num2str(handles.gonad.im4_scale));
    
    set(handles.im1_off,  'String',num2str(handles.gonad.im1_off));
    set(handles.im2_off,  'String',num2str(handles.gonad.im2_off));
    set(handles.im3_off,  'String',num2str(handles.gonad.im3_off));
    set(handles.im4_off,  'String',num2str(handles.gonad.im4_off));
    
    set(handles.im1_gamma,'String',num2str(handles.gonad.im1_gamma));
    set(handles.im2_gamma,'String',num2str(handles.gonad.im2_gamma));
    set(handles.im3_gamma,'String',num2str(handles.gonad.im3_gamma));
    set(handles.im4_gamma,'String',num2str(handles.gonad.im4_gamma));
    
end

if isfield(handles.gonad, 'im1_index')
    handles.im1_index = handles.gonad.im1_index;
    handles.im2_index = handles.gonad.im2_index;
    handles.im3_index = handles.gonad.im3_index;
    handles.im4_index = handles.gonad.im4_index;
else
    handles.im1_index = 1;
    handles.im2_index = 2;
    handles.im3_index = 3;
    handles.im4_index = [];
end

if isempty(intersect(handles.im1_index, 1))
    set(handles.button_im1_R, 'FontWeight', 'Normal');
else
    set(handles.button_im1_R, 'FontWeight', 'Bold');
end
if isempty(intersect(handles.im1_index, 2))
    set(handles.button_im1_G, 'FontWeight', 'Normal');
else
    set(handles.button_im1_G, 'FontWeight', 'Bold');
end
if isempty(intersect(handles.im1_index, 3))
    set(handles.button_im1_B, 'FontWeight', 'Normal');
else
    set(handles.button_im1_B, 'FontWeight', 'Bold');
end


if isempty(intersect(handles.im2_index, 1))
    set(handles.button_im2_R, 'FontWeight', 'Normal');
else
    set(handles.button_im2_R, 'FontWeight', 'Bold');
end
if isempty(intersect(handles.im2_index, 2))
    set(handles.button_im2_G, 'FontWeight', 'Normal');
else
    set(handles.button_im2_G, 'FontWeight', 'Bold');
end
if isempty(intersect(handles.im2_index, 3))
    set(handles.button_im2_B, 'FontWeight', 'Normal');
else
    set(handles.button_im2_B, 'FontWeight', 'Bold');
end


if isempty(intersect(handles.im3_index, 1))
    set(handles.button_im3_R, 'FontWeight', 'Normal');
else
    set(handles.button_im3_R, 'FontWeight', 'Bold');
end
if isempty(intersect(handles.im3_index, 2))
    set(handles.button_im3_G, 'FontWeight', 'Normal');
else
    set(handles.button_im3_G, 'FontWeight', 'Bold');
end
if isempty(intersect(handles.im3_index, 3))
    set(handles.button_im3_B, 'FontWeight', 'Normal');
else
    set(handles.button_im3_B, 'FontWeight', 'Bold');
end


if isfield(handles.gonad, 'chan1_tag')
    set(handles.chan1_tag, 'string', handles.gonad.chan1_tag);
    set(handles.chan2_tag, 'string', handles.gonad.chan2_tag);
    set(handles.chan3_tag, 'string', handles.gonad.chan3_tag);
    set(handles.chan4_tag, 'string', handles.gonad.chan4_tag);
end

try
    set(handles.seg_thresh,         'String', num2str(handles.gonad.SEG_THRESH));
    set(handles.seg_channel,        'String', num2str(handles.gonad.SEG_CHANNEL));
    set(handles.seg_smoothing,      'String', num2str(handles.gonad.SEG_SMOOTHING));
    set(handles.seg_gamma,          'String', num2str(handles.gonad.SEG_GAMMA));
end

try
    set(handles.checkbox_useSegGamma,   'Value', handles.gonad.useGamma);
    set(handles.checkbox_useOtsuThresh, 'Value', handles.gonad.useOtsuThresh);
end

handles.focus_channel_1_status = 1;
handles.focus_channel_2_status = 0;
set(handles.button_focus_channel_1, 'FontWeight', 'Bold');

% load default focus-axis correspondence
if isfield(handles.gonad, 'focusAxisTable')
    dat = handles.gonad.focusAxisTable;
else
    dat = zeros(6, 2);
end

set(handles.focusAxisTable, 'Data', dat);
set(handles.focusAxisTable, 'ColumnEditable', logical(ones(1, size(dat,2))));

if isfield(handles.gonad, 'focusChannels')
    set(handles.focus_channel_1,    'String', num2str(handles.gonad.focusChannels(1)));
    set(handles.focus_channel_2,    'String', num2str(handles.gonad.focusChannels(2)));
end

if isfield(handles.gonad, 'focusSmoothing')
    set(handles.focus_smoothing_1,  'String', num2str(handles.gonad.focusSmoothing(1)));
    set(handles.focus_smoothing_2,  'String', num2str(handles.gonad.focusSmoothing(2)));
end

if isfield(handles.gonad, 'focusThresh')
    set(handles.focus_thresh_1,     'String', num2str(handles.gonad.focusThresh(1)));
    set(handles.focus_thresh_2,     'String', num2str(handles.gonad.focusThresh(2)));
end

handles.nukeDir = 'nukes';

nukeFolders = dir([handles.dd filesep handles.nukeDir ]);
counter = 0;
good = [];

for ii = 1:length(nukeFolders)
    
    if (nukeFolders(ii).isdir)
        
        if  isempty(intersect(str2double(nukeFolders(ii).name), 1:999))
            continue;
        end
        
        counter = counter +1;
        nukeFoldersList{counter} = nukeFolders(ii).name;
        
        good = [good, str2double(nukeFolders(ii).name)];
        
    end
    
end

handles.nukeFoldersList = nukeFoldersList;
handles.gonad.good_nukes = good;

if ~isfield(handles.gonad, 'nuke_ids') || isempty(handles.gonad.nuke_ids)
    handles.gonad.nuke_ids = good*0 + 3;
end

if ~isfield(handles.gonad, 'finalized_nukes') || isempty(handles.gonad.finalized_nukes)
    handles.gonad.finalized_nukes = good*0;
end

set(handles.textbox_number_finalize, 'string', num2str(sum(handles.gonad.finalized_nukes==1)));
set(handles.textbox_number_error,    'string', num2str(sum(handles.gonad.finalized_nukes==2)));
set(handles.textbox_number_ignore,   'string', num2str(sum(handles.gonad.finalized_nukes==-1)));


handles.index_selected = 1;
set(handles.listbox_surfs, 'String', handles.nukeFoldersList, 'Value', 1);

% Nucleus labels and projection type selection (mask or image)
handles.showLabels              = 0;
handles.useGonadMaskProjection  = 0;

set(handles.checkbox_use_smoothing, 'Value', 0);

if handles.useGonadMaskProjection
    [handles.mkF, handles.mkFC_projection] = load_gonad_mask_projection(handles);
    handles.mkFC_projection_orig = handles.mkFC_projection;
else
    handles.imProjection = load_gonad_image_projection(handles);
end

if ~isfield(handles.gonad, 'nukeProps')

    if ~handles.useGonadMaskProjection
        [handles.mkF, handles.mkFC_projection] = load_gonad_mask_projection(handles);
    end
    
    handles.gonad.nukeProps = regionprops(handles.mkF, 'Centroid','BoundingBox');

end

if handles.useGonadMaskProjection

    for ii = 1:length(handles.gonad.finalized_nukes)

        % ignored nucleus
        if handles.gonad.finalized_nukes(ii)==-1
            cc = [155,155,155];

        % good nucleus
        elseif handles.gonad.finalized_nukes(ii)==1
            cc = [55,255,55];

        % error nucleus
        elseif handles.gonad.finalized_nukes(ii)==2
            cc = [255,55,0];

        end

        imb = handles.mkFC_projection;
        mk = max(handles.mkF==handles.gonad.good_nukes(ii),[],3);

        imR = imb(:,:,1);   imR(mk) = cc(1);    handles.mkFC_projection(:,:,1) = imR;
        imG = imb(:,:,2);   imG(mk) = cc(2);    handles.mkFC_projection(:,:,2) = imG;
        imB = imb(:,:,3);   imB(mk) = cc(3);    handles.mkFC_projection(:,:,3) = imB;

    end

end


handles.labelAxesInVolumeOrder = 0;

if isfield(handles.gonad, 'labelAxesInVolumeOrder')
    handles.labelAxesInVolumeOrder = handles.gonad.labelAxesInVolumeOrder;
end

end

% ============================================================--
%
%  PROJECTED MASK LOADING
%
% ============================================================--


function [mkF, mkFC_projection] = load_gonad_mask_projection(handles)

mkF = load3([handles.dd filesep 'mkF3_reseg.tif']);

try
    mkFC_ = imread([handles.dd filesep 'mkF3_reseg_color_proj.tif']);
catch
    mkFC = makeColorMask3(double(mkF));
    mkFC_ = squeeze(max(mkFC,[],3));
    imwrite(mkFC_, [handles.dd filesep 'mkF3_reseg_color_proj.tif']);
end

tmp = max(mkFC_,[],3)==0;
mkFC_ = mkFC_ + uint8(0);

mkFC_(cat(3, tmp,tmp,tmp)) = 0;
mkFC_projection = mkFC_;

end

function imProjection = load_gonad_image_projection(handles)


for ii = 1:4
    
    if ii <= handles.gonad.NUM_CHAN
        
        eval(['im{' num2str(ii) '} = imread([handles.dd filesep ''im' num2str(ii) '_PROJ.tif'']);']);

        
        if size(size(im{ii}),2)>1
            im{ii} = im{ii}(:,:,1);
        end
        
    else
        
        im{ii} = 0*im{1};
        
    end
end
    
scale = [str2double(get(handles.im1_scale,'String')),...
         str2double(get(handles.im2_scale,'String')),...
         str2double(get(handles.im3_scale,'String')),...
         str2double(get(handles.im4_scale,'String'))];

offset = [str2double(get(handles.im1_off,'String')),...
          str2double(get(handles.im2_off,'String')),...
          str2double(get(handles.im3_off,'String')),...
          str2double(get(handles.im4_off,'String'))];

gamma = [str2double(get(handles.im1_gamma,'String')),...
         str2double(get(handles.im2_gamma,'String')),...
         str2double(get(handles.im3_gamma,'String')),...
         str2double(get(handles.im4_gamma,'String'))];

     scale = scale;
     % gamma = [1,1,1,1];

for ii = 1:length(im)

    im{ii}  = (scale(ii)*(im{ii} - mean(im{ii}(:))*offset(ii)));
    im{ii}  = autogain(double(im{ii}).^(gamma(ii)));

end


imr = autogain(0*im{1});
img = autogain(0*im{1});
imb = autogain(0*im{1});


if sum(handles.im1_index==1)
    imr = imr + im{1};
end
if sum(handles.im1_index==2)
    img = img + im{1};
end
if sum(handles.im1_index==3)
    imb = imb + im{1};
end
if sum(handles.im2_index==1)
    imr = imr + im{2};
end
if sum(handles.im2_index==2)
    img = img + im{2};
end
if sum(handles.im2_index==3)
    imb = imb + im{2};
end
if sum(handles.im3_index==1)
    imr = imr + im{3};
end
if sum(handles.im3_index==2)
    img = img + im{3};
end
if sum(handles.im3_index==3)
    imb = imb + im{3};
end
if sum(handles.im4_index==1)
    imr = imr + im{4};
end
if sum(handles.im4_index==2)
    img = img + im{4};
end
if sum(handles.im4_index==3)
    imb = imb + im{4};
end


rr = [200, 25, 25]/255;
gg = [25, 200, 25]/255;
bb = [25, 25, 200]/255;


imProjection = cat(3,...
    imr*rr(1) + img*gg(1) + imb*bb(1),...
    imr*rr(2) + img*gg(2) + imb*bb(2),...
    imr*rr(3) + img*gg(3) + imb*bb(3));




end


% --- Executes on button press in vizIm.
function vizIm_Callback(hObject, eventdata, handles)

display_nucleus_vol_viz(handles);
try
    display_axis_images(handles);
end

end

% ============================================================--
%
%  NUCLEUS IDENTITY BUTTONS (ARCHAIC)
%
% ============================================================--


% --- Executes on button press in button_pachy.
function button_pachy_Callback(hObject, eventdata, handles)

if handles.gonad.nuke_ids(handles.index_selected) == 3
    set(handles.button_pachy, 'FontWeight', 'Normal');
    handles.gonad.nuke_ids(handles.index_selected) = 0;
else
    set(handles.button_tz, 'FontWeight', 'Normal');
    set(handles.button_earlyPachy, 'FontWeight', 'Normal');
    set(handles.button_pachy, 'FontWeight', 'Bold');
    handles.gonad.nuke_ids(handles.index_selected) = 3;
end
display_gonad(handles, hObject);
guidata(hObject, handles);

end


% --- Executes on button press in button_earlyPachy.
function button_earlyPachy_Callback(hObject, eventdata, handles)

if handles.gonad.nuke_ids(handles.index_selected) == 2
    set(handles.button_earlyPachy, 'FontWeight', 'Normal');
    handles.gonad.nuke_ids(handles.index_selected) = 0;
else
    set(handles.button_tz, 'FontWeight', 'Normal');
    set(handles.button_earlyPachy, 'FontWeight', 'Bold');
    set(handles.button_pachy, 'FontWeight', 'Normal');
    handles.gonad.nuke_ids(handles.index_selected) = 2;
end
display_gonad(handles, hObject);
guidata(hObject, handles);

end

% --- Executes on button press in button_tz.
function button_tz_Callback(hObject, eventdata, handles)

if handles.gonad.nuke_ids(handles.index_selected) == 1
    set(handles.button_tz, 'FontWeight', 'Normal');
    handles.gonad.nuke_ids(handles.index_selected) = 0;
else
    set(handles.button_tz, 'FontWeight', 'Bold');
    set(handles.button_earlyPachy, 'FontWeight', 'Normal');
    set(handles.button_pachy, 'FontWeight', 'Normal');
    handles.gonad.nuke_ids(handles.index_selected) = 1;
end
display_gonad(handles, hObject);
guidata(hObject, handles);
end

% --- Executes on button press in button_clearIDs.
function button_clearIDs_Callback(hObject, eventdata, handles)

handles.gonad.nuke_ids = zeros(length(handles.nukeFoldersList),1);

set(handles.button_tz, 'FontWeight', 'Normal');
set(handles.button_earlyPachy, 'FontWeight', 'Normal');
set(handles.button_pachy, 'FontWeight', 'Normal');

display_gonad(handles, hObject);
guidata(hObject, handles);

end


% ============================================================--
%
%  NUCLEUS 3D IMAGE DISPLAY
%
% ============================================================--


function display_nucleus_vol_viz(handles, figNum, use_resz)

if ~exist('figNum', 'var')
    figNum = 2;
end

if ~exist('use_resz', 'var')
    use_resz = 0;
end

figure(2); clf;

[im, channels] = make_image(handles, use_resz);

if sum(~~channels)==1 && get(handles.checkbox_use_colormap, 'Value')
    
    im = squeeze(im(:,:,:,channels(find(~~channels))));
    
    vizIm(autogain((double(im)+1)), figNum);
    
    cc = colormap(jet(255));
    colormap(flipud(cc));
    
else
    vizIm(im, figNum);
end

ff = figure(figNum);
handles.alphamap = str2double(get(handles.alpha, 'String'));
set(ff, 'Alphamap', (0:64)*handles.alphamap/64);

screen_sz = getScreenSize;

if screen_sz(3)==2560
    
    set(ff,'position',  [780,  100,  1050,  750]);
   
% if home monitor
elseif screen_sz(3)==1920
    
    set(ff, 'position', [750,  580,  580,   370]);
    
end

rotate3d;


% set(ff, 'WindowKeyPressFcn', @presser);

    function presser(src, event)
        
        if char(event.Character) == 'i'
            button_ignore_Callback(hObject, [], handles);
        end
        
        if char(event.Character) == 'f'
            button_finalize_Callback(hObject, [], handles)
        end
        
        if char(event.Character) == 'e'
            button_error_Callback(hObject, [], handles)
        end
        
    end

end

% ============================================================--
%
%  NUCLEUS IMAGE CONSTRUCTION (INTERNAL)
%
% ============================================================--


function [im, channels] = make_image(handles, use_resz)

try
    mk = uint8(imdilate(handles.mk, ones(11,11,11)));
end

scale = [str2double(get(handles.im1_scale,'String')),...
         str2double(get(handles.im2_scale,'String')),...
         str2double(get(handles.im3_scale,'String')),...
         str2double(get(handles.im4_scale,'String'))];

offset = [str2double(get(handles.im1_off,'String')),...
          str2double(get(handles.im2_off,'String')),...
          str2double(get(handles.im3_off,'String')),...
          str2double(get(handles.im4_off,'String'))];

gamma = [str2double(get(handles.im1_gamma,'String')),...
         str2double(get(handles.im2_gamma,'String')),...
         str2double(get(handles.im3_gamma,'String')),...
         str2double(get(handles.im4_gamma,'String'))];
     
if 1 % use_resz
    im1 = (scale(1)*(handles.im1 - mean(handles.im1(:))*offset(1)));
    im2 = (scale(2)*(handles.im2 - mean(handles.im2(:))*offset(2)));
    im3 = (scale(3)*(handles.im3 - mean(handles.im3(:))*offset(3)));
    im4 = (scale(4)*(handles.im4 - mean(handles.im4(:))*offset(4)));
else
    im1 = (scale(1)*(handles.im1_raw - mean(handles.im1_raw(:))*offset(1)));
    im2 = (scale(2)*(handles.im2_raw - mean(handles.im2_raw(:))*offset(2)));
    im3 = (scale(3)*(handles.im3_raw - mean(handles.im3_raw(:))*offset(3)));
    im4 = (scale(4)*(handles.im4_raw - mean(handles.im4_raw(:))*offset(4)));
end

im1  = autogain(double(im1).^(gamma(1)));
im2  = autogain(double(im2).^(gamma(2)));
im3  = autogain(double(im3).^(gamma(3)));
im4  = autogain(double(im4).^(gamma(4)));

ord = [handles.im1_index, handles.im2_index, handles.im3_index, handles.im4_index];

channels = (ord);

r_inds = find(ord==1);
g_inds = find(ord==2);
b_inds = find(ord==3);

imr = autogain(0*im1);
img = autogain(0*im1);
imb = autogain(0*im1);

% for ii = r_inds
%     imr = imr + eval(['im' num2str(ii) ';']);
% end
% for ii = g_inds
%     img = img + eval(['im' num2str(ii) ';']);
% end
% for ii = b_inds
%     imb = imb + eval(['im' num2str(ii) ';']);
% end

if sum(handles.im1_index==1)
    imr = imr + im1;
end
if sum(handles.im1_index==2)
    img = img + im1;
end
if sum(handles.im1_index==3)
    imb = imb + im1;
end
if sum(handles.im2_index==1)
    imr = imr + im2;
end
if sum(handles.im2_index==2)
    img = img + im2;
end
if sum(handles.im2_index==3)
    imb = imb + im2;
end
if sum(handles.im3_index==1)
    imr = imr + im3;
end
if sum(handles.im3_index==2)
    img = img + im3;
end
if sum(handles.im3_index==3)
    imb = imb + im3;
end
if sum(handles.im4_index==1)
    imr = imr + im4;
end
if sum(handles.im4_index==2)
    img = img + im4;
end
if sum(handles.im4_index==3)
    imb = imb + im4;
end

if get(handles.checkbox_useMask, 'Value')
    
    try
        imr = imr.*mk;
        img = img.*mk;
        imb = imb.*mk;
    catch
        'WARNING: nucleus masking not working'
    end
    
    
end

rr = [255, 0, 60]/255;
gg = [60, 255, 0]/255;
bb = [0, 125, 255]/255;

rr = [200, 25, 25]/255;
gg = [25, 200, 25]/255;
bb = [25, 25, 200]/255;

% rr = [1, 1, 0];
% gg = [1, 0, 1];
% bb = [0, 1, 1];
%
% rr = [1, 1/2, 1/2];
% gg = [1/2, 1, 1/2];
% bb = [1/2, 1/2, 1];

imR = imr*rr(1) + img*gg(1) + imb*bb(1);
imG =  imr*rr(2) + img*gg(2) + imb*bb(2);
imB = imr*rr(3) + img*gg(3) + imb*bb(3);

im = cat(4,...
    imr*rr(1) + img*gg(1) + imb*bb(1),...
    imr*rr(2) + img*gg(2) + imb*bb(2),...
    imr*rr(3) + img*gg(3) + imb*bb(3));



% im = cat(4,autogain(imr),autogain(img),autogain(imb));

end

% ============================================================--
%
%  GONAD DISPLAY USING MASK
%
% ============================================================--

function display_gonad(handles, hObject)

ff = figure(1);clf;
set(ff, 'WindowButtonDownFcn', @clicker);
set(ff, 'WindowKeyPressFcn', @presser);

set(gca,'Position',[0,0,1,1]);
screen_sz = getScreenSize;

% if 30 inch monitor
if screen_sz(3)==2560
    if handles.gonad.imageSize(1) > 512
        set(gcf, 'position', [1550           55        1000    1000]);
    else
        set(gcf, 'position', [2560-770-100   1600-850   850     770]);
    end
    
end

% if home monitor
if screen_sz(3)==1920
    
    set(gcf,'position',      [1050            700       850      770]);
    
end


display_image_internal(handles);
guidata(hObject, handles);

    function display_image_internal(handles)
        
        if handles.useGonadMaskProjection
    
            mk = max(handles.mkF==handles.gonad.good_nukes(handles.index_selected),[],3);
        
            mk = ~~(imdilate(~~mk, ones(4,4,4)) - (~~mk));

            imb = handles.mkFC_projection;

            imR = imb(:,:,1);
            imG = imb(:,:,2);
            imB = imb(:,:,3);

            imR(mk) = 255;
            imG(mk) = 255;
            imB(mk) = 255;


            imshow(cat(3, imR, imG, imB),[]);
            hold on
        
        else
            
            imb = handles.imProjection;
            
            imshow(imb,[]);
            hold on
            
            xx = cos(0:.1:6.5);
            yy = sin(0:.1:6.5);
                        
            for ii = 1:length(handles.gonad.good_nukes)
                
                [rr,cc,zz] = getBoundingBox3(handles.gonad.nukeProps, handles.gonad.good_nukes(ii), 0, handles.gonad.imageSize);
                
                pos = handles.gonad.nukeProps(handles.gonad.good_nukes(ii)).Centroid;
                
                rad = [max(rr) - min(rr), max(cc) - min(cc)]/2;
                
                % ignored nucleus
                if handles.gonad.finalized_nukes(ii)==-1
                    circle_color = [1,1,1]/3;

                % good nucleus
                elseif handles.gonad.finalized_nukes(ii)==1
                    circle_color = [55,255,55]/255;

                % error nucleus
                elseif handles.gonad.finalized_nukes(ii)==2
                    circle_color = [255,55,0]/255;

                % no identity
                else
                    circle_color = [55, 55, 255]/255;
                end
                
                lineWidth = 2;
                
                if ii==handles.index_selected
                    lineWidth = 4;
                    circle_color = circle_color + [1,1,1]/2;
                    circle_color = circle_color / max(circle_color);
                end
                
                plot(xx*rad(2) + pos(1), yy*rad(1) + pos(2), 'color', circle_color, 'linewidth', lineWidth);

            end   
        end   
        
        % this runs whether we're using the mask or image projection
        if handles.showLabels
            
            for ii = 1:length(handles.gonad.good_nukes)
                
                % transition zone
                if handles.gonad.nuke_ids(ii)==1
                    cc = [144, 238, 144]/255;

                % early pachytene
                elseif handles.gonad.nuke_ids(ii)==2
                    cc = [130,200,240]/255;
                                        
                % pachytene
                elseif handles.gonad.nuke_ids(ii)==3
                    cc = [70,100,240]/255;
                        
                % unlabeled
                else
                    cc = [1,1,1];
                end


                text(...
                    handles.gonad.nukeProps(handles.gonad.good_nukes(ii)).Centroid(1),...
                    handles.gonad.nukeProps(handles.gonad.good_nukes(ii)).Centroid(2),...
                    num2str(handles.gonad.good_nukes(ii)), 'fontsize', 36,'fontweight', 'bold', 'color', cc);
                
            end
        end
    end


% ============================================================--
%
%  ONCLICK FUNCTION FOR NUCLEUS SELECTION
%
% ============================================================--

    function clicker(src, event)
        
        ord = [handles.im1_index, handles.im2_index, handles.im3_index];
        
        click_point = get(gca,'CurrentPoint');
        click_point = click_point(1,1:2);
        
        for ii = 1:length(handles.gonad.good_nukes)
            
            dist(ii) = sum((handles.gonad.nukeProps(handles.gonad.good_nukes(ii)).Centroid(1:2) - click_point).^2);
            
        end
        
        minn = min(dist);
        ind = find(dist==minn);
        
        if ind == handles.index_selected
            return;
        end
        
        handles.index_selected = ind;
        handles.gonad.good_nukes(ind)
        
        nd = load([
            handles.dd filesep ...
            handles.nukeDir filesep ...
            num2str(handles.gonad.good_nukes(ind)) filesep ...
            'nd.mat']);
        
        nd = nd.nd;
        
        nd.PX_SZ = handles.gonad.PHYSICAL_PX_SZ(1);
        
        scale = [str2double(get(handles.im1_scale,'String')),...
                 str2double(get(handles.im2_scale,'String')),...
                 str2double(get(handles.im3_scale,'String')),...
                 str2double(get(handles.im4_scale,'String'))]; 
        
        offset = [str2double(get(handles.im1_off,'String')),...
                  str2double(get(handles.im2_off,'String')),...
                  str2double(get(handles.im3_off,'String')),...
                  str2double(get(handles.im4_off,'String'))];
        
        gamma = [str2double(get(handles.im1_gamma,'String')),...
                 str2double(get(handles.im2_gamma,'String')),...
                 str2double(get(handles.im3_gamma,'String')),...
                 str2double(get(handles.im4_gamma,'String'))];
        
        handles.nd.displaySettings.scale   = scale;
        handles.nd.displaySettings.offset  = offset;
        handles.nd.displaySettings.gamma   = gamma;
        
        handles.nd = nd;
        handles.resave_nd_flag = 0;
        
        handles.im1      = handles.nd.im1;
        handles.im1_raw  = handles.nd.im1;
        
        if handles.gonad.NUM_CHAN>=2
            handles.im2     = handles.nd.im2;
            handles.im2_raw = handles.nd.im2;
        else
            handles.im2     = [];
            handles.im2_raw = [];
        end
        
        if handles.gonad.NUM_CHAN>=3
            handles.im3     = handles.nd.im3;
            handles.im3_raw = handles.nd.im3;
        else
            handles.im3     = [];
            handles.im3_raw = [];
        end
        
        if handles.gonad.NUM_CHAN==4
            handles.im4     = handles.nd.im4;
            handles.im4_raw = handles.nd.im4;
        else
            handles.im4     = [];
            handles.im4_raw = [];
        end
        
        handles.mk = handles.nd.mk;

        
        if get(handles.checkbox_use_smoothing,'Value')
            
            RESZ = 2;
            
            handles.im1         = imresize3(handles.nd.im1, RESZ);
            handles.mk          = imresize3(handles.mk, RESZ);
            
            if ~isempty(handles.im2)
                handles.im2     = imresize3(handles.nd.im2, RESZ);
            else
                handles.im2     = handles.im1;
                handles.im2_raw = handles.im1_raw;
            end
            
            if ~isempty(handles.im3)
                handles.im3     = imresize3(handles.nd.im3, RESZ);
            else
                handles.im3     = handles.im1;
                handles.im3_raw = handles.im1_raw;
            end
            
            if ~isempty(handles.im4)
                handles.im4     = imresize3(handles.nd.im4, RESZ);
            else
                handles.im4     = handles.im1;
                handles.im4_raw = handles.im1_raw;
            end
            
        else
            
            if isempty(handles.im2)
                handles.im2 = handles.im1;
            end
            
            if isempty(handles.im3)
                handles.im3 = handles.im1;
            end
            
            if isempty(handles.im4)
                handles.im4 = handles.im1;
            end
            
        end
        
        if isfield(nd, 'segParams')
            if isfield(nd.segParams, 'GAMMA')
                nd.segParams.SEG_GAMMA = nd.segParams.GAMMA;
            end
            
            try
                set(handles.seg_thresh,     'String', num2str(nd.segParams.SEG_THRESH));
                set(handles.seg_channel,    'String', num2str(nd.segParams.SEG_CHANNEL));
                set(handles.seg_smoothing,  'String', num2str(nd.segParams.SEG_SMOOTHING));
                set(handles.seg_gamma,      'String', num2str(nd.segParams.SEG_GAMMA));
            catch
                'WARNING: failed to load seg settings'
            end
            
            try
                set(handles.checkbox_useOtsuThresh, 'Value', nd.segParams.useOtsuThresh);
                set(handles.checkbox_useSegGamma, 'Value', nd.segParams.useGamma);
            catch
                'WARNING: failed to load gamma/thresh settings'
            end
            
        end
        
        if isfield(nd, 'focusFittingParams')
            
            len = length(nd.focusFittingParams.channel);
            
            if len==2
                try
                    set(handles.focus_smoothing_1, 'String', num2str(nd.focusFittingParams.smoothing(1)));
                    set(handles.focus_smoothing_2, 'String', num2str(nd.focusFittingParams.smoothing(2)));
                end
                try
                    set(handles.focus_thresh_1,  'String', num2str(nd.focusFittingParams.thresh(1)));
                    set(handles.focus_thresh_2,  'String', num2str(nd.focusFittingParams.thresh(2)));
                end
            else
                try
                    set(handles.focus_smoothing_1, 'String', num2str(nd.focusFittingParams.smoothing(1)));
                end
                try
                    set(handles.focus_thresh_1,  'String', num2str(nd.focusFittingParams.thresh(1)));
                end
            end
        end
        
        
        % this channel assignment takes precedence if present
        % since zeros here prevent foci from that channel from being used by optiAxis
        if isfield(nd, 'focusChannel')
            
            set(handles.focus_channel_1, 'String', num2str(nd.focusChannel(1)));
            
            if length(nd.focusChannel)>1
                set(handles.focus_channel_2, 'String', num2str(nd.focusChannel(2)));
            end
            
        elseif isfield(nd, 'focusFittingParams')
            
            set(handles.focus_channel_1, 'String', num2str(nd.focusFittingParams.channel(1)));
            
            if length(nd.focusFittingParams.channel)>1
                set(handles.focus_channel_2, 'String', num2str(nd.focusFittingParams.channel(2)));
            end
            
        end
        
        
        
        if isfield(nd, 'focusAxisTable')
            set(handles.focusAxisTable, 'data', nd.focusAxisTable);
        else
            % load default focus-axis correspondence
            if isfield(handles.gonad, 'focusAxisTable')
                dat = handles.gonad.focusAxisTable;
            else
                dat = zeros(6, 2);
            end
            
            set(handles.focusAxisTable, 'Data', dat);
        end
        
        
        try
            if strcmp(nd.optParams.OPTZ_METHOD, 'intensity')
                set(handles.radio_intensity, 'Value', get(handles.radio_intensity,'Max'));
            elseif strcmp(nd.optParams.OPTZ_METHOD, 'random')
                set(handles.radio_random, 'Value', get(handles.radio_random,'Max'));
            elseif strcmp(nd.optParams.OPTZ_METHOD, 'orientation')
                set(handles.radio_orientation, 'Value', get(handles.radio_orientation,'Max'));
            end
        end
        
        
        set(handles.button_tz,          'FontWeight', 'Normal');
        set(handles.button_earlyPachy,  'FontWeight', 'Normal');
        set(handles.button_pachy,       'FontWeight', 'Normal');
        
        set(handles.button_ignore,      'FontWeight', 'Normal');
        set(handles.button_finalize,    'FontWeight', 'Normal');
        set(handles.button_error,       'FontWeight', 'Normal');
        
        if handles.gonad.nuke_ids(handles.index_selected) == 3
            set(handles.button_pachy, 'FontWeight', 'Bold');
        elseif handles.gonad.nuke_ids(handles.index_selected) == 2
            set(handles.button_earlyPachy, 'FontWeight', 'Bold');
        elseif handles.gonad.nuke_ids(handles.index_selected) == 1
            set(handles.button_tz, 'FontWeight', 'Bold');
        end
        
        if handles.gonad.finalized_nukes(handles.index_selected) == 1
            set(handles.button_finalize, 'FontWeight', 'Bold');
        elseif handles.gonad.finalized_nukes(handles.index_selected) == -1
            set(handles.button_ignore, 'FontWeight', 'Bold');
        elseif handles.gonad.finalized_nukes(handles.index_selected) == 2
            set(handles.button_error, 'FontWeight', 'Bold');
        end
        
        
        display_image_internal(handles);
        display_nucleus_vol_viz(handles);
        
        status = [];
        
        if isfield(handles.nd, 'optiData')
            
            status = handles.nd.optiData.status;
            
        elseif isfield(handles.nd, 'surf_status')
            
            status = handles.nd.surf_status;
            
        else
            if isfield(handles.nd, 'sdata')
                nd.optiData.status = zeros(1, length(nd.sdata));
            end
        end
        
        if ~isempty(status)
            
            [on,off] = make_axis_mask(nd.surfL, status);
            
            if get(handles.checkbox_displaySurfs,'value')
                show_labels = get(handles.checkbox_displaySurfLabels, 'value');
                vizSegInitial(handles.nd, on, off, 11, 1, show_labels);
            end
            
            if get(handles.checkbox_displayFinalSeg, 'value')
                
                handles.nd = makeAxisMask(handles.nd, handles.labelAxesInVolumeOrder);
                
                figNum = 111;
                handles.nd = displayFinalSeg(handles.nd, figNum);
                
                screen_sz = getScreenSize;
                
                % if 30 inch monitor
                if screen_sz(3)==2560
                    set(gcf, 'position', [850, 1050, 560, 420]);
                end

                handles = load_axis_mask_status(handles);
                
            end
            
            % add all possible surfs to the list box
            handles.inds = 1:length(handles.nd.ppS);
            vals = [];
            
            for ii = 1:length(handles.inds)
                inds_str{ii} = num2str(handles.inds(ii));
                if status(handles.inds(ii))
                    vals = [vals, ii];
                end
            end
            
            set(handles.listbox_surfs, 'String', inds_str, 'Value', vals);
            
            handles.vals = vals;
            handles.inds_str = inds_str;
            
        end
        
        guidata(hObject, handles);
        figure(1);
        
    end

    function presser(src, event)
        
        if char(event.Character) == 'd'
            handles.do_split = 0;
            button_ignore_Callback(hObject, [], handles);
        end
        
        if char(event.Character) == 'f'
            button_finalize_Callback(hObject, [], handles)
        end
        
        if char(event.Character) == 'e'
            button_error_Callback(hObject, [], handles)
        end
        
    end


figure(1);
end




% ============================================================--
%
%  SAVE SETTINGS BUTTON
%
% ============================================================--


% --- Executes on button press in button_saveSettings.
function button_saveSettings_Callback(hObject, eventdata, handles)

handles.gonad.im1_scale = get(handles.im1_scale,'String');
handles.gonad.im2_scale = get(handles.im2_scale,'String');
handles.gonad.im3_scale = get(handles.im3_scale,'String');
handles.gonad.im4_scale = get(handles.im4_scale,'String');

handles.gonad.im1_off = get(handles.im1_off,'String');
handles.gonad.im2_off = get(handles.im2_off,'String');
handles.gonad.im3_off = get(handles.im3_off,'String');
handles.gonad.im4_off = get(handles.im4_off,'String');


handles.gonad.im1_gamma = get(handles.im1_gamma,'String');
handles.gonad.im2_gamma = get(handles.im2_gamma,'String');
handles.gonad.im3_gamma = get(handles.im3_gamma,'String');
handles.gonad.im4_gamma = get(handles.im4_gamma,'String');

handles.gonad.im1_index = handles.im1_index;
handles.gonad.im2_index = handles.im2_index;
handles.gonad.im3_index = handles.im3_index;
handles.gonad.im4_index = handles.im4_index;

gonad = handles.gonad;

gonad.SEG_THRESH     = str2double(get(handles.seg_thresh, 'String'));
gonad.SEG_CHANNEL    = str2double(get(handles.seg_channel, 'String'));
gonad.SEG_SMOOTHING  = str2double(get(handles.seg_smoothing, 'String'));
gonad.SEG_GAMMA      = str2double(get(handles.seg_gamma, 'String'));

gonad.useGamma = get(handles.checkbox_useSegGamma, 'Value');
gonad.useOtsuThresh = get(handles.checkbox_useOtsuThresh, 'Value');

if (get(handles.radio_intensity,'Value') == get(handles.radio_intensity,'Max'))
    order_param_str = 'intensity';
end
if (get(handles.radio_orientation,'Value') == get(handles.radio_orientation,'Max'))
    order_param_str = 'orientation';
end
if (get(handles.radio_random,'Value') == get(handles.radio_random,'Max'))
    order_param_str = 'random';
end

% gonad.MAX_AXIS_LENGTH = str2double(get(handles.seg_gamma, 'String'));
%
% gonad.OPTZ_METHOD = order_param_str;

save( [handles.dd filesep 'gonad.mat'], 'gonad');

[handles.dd filesep 'gonad.mat']

end

% ============================================================--
%
%  SHOW THRESHOLD BUTTON
%
% ============================================================--


% --- Executes on button press in button_showThresh.
function button_showThresh_Callback(hObject, eventdata, handles)

SEG_CHANNEL      = str2double(get(handles.seg_channel, 'String'));
SEG_THRESH       = str2double(get(handles.seg_thresh, 'String'));
SEG_SMOOTHING    = str2double(get(handles.seg_smoothing, 'String'));
SEG_GAMMA        = str2double(get(handles.seg_gamma, 'String'));

% set these values to pass to makeAxisBackMask, but do not save them
% so that they do not overwrite existing values unless the segmentation
% itself is re-run.

handles.nd.segParams.SEG_CHANNEL     = SEG_CHANNEL;
handles.nd.segParams.SEG_THRESH      = SEG_THRESH;
handles.nd.segParams.SEG_SMOOTHING   = SEG_SMOOTHING;
handles.nd.segParams.SEG_GAMMA       = SEG_GAMMA;

handles.nd.segParams.useGamma       = get(handles.checkbox_useSegGamma, 'Value');
handles.nd.segParams.useOtsuThresh  = get(handles.checkbox_useOtsuThresh, 'Value');

handles.nd = makeAxisBackMask(handles.nd);

% mkAxisBack = mkAxisBack .* imdilate(mkNuke,ones(5,5,5));

vizMask(handles.nd.mkAxis, [1,1,1],8,1,1,1);

set(handles.seg_thresh, 'String', num2str(handles.nd.segParams.SEG_THRESH));


end

% ============================================================--
%
%  SEGMENTATION BUTTON
%
% ============================================================--


% --- Executes on button press in button_segment.
function button_segment_Callback(hObject, eventdata, handles)

'==============================  SEGMENTING =============================='


handles.resave_nd_flag = 1;

SEG_CHANNEL      = str2double(get(handles.seg_channel, 'String'));
SEG_THRESH       = str2double(get(handles.seg_thresh, 'String'));
SEG_SMOOTHING    = str2double(get(handles.seg_smoothing, 'String'));
SEG_GAMMA        = str2double(get(handles.seg_gamma, 'String'));

% set these values to pass to makeAxisBackMask, but do not save them
% so that they do not overwrite existing values unless the segmentation
% itself is re-run.

handles.nd.segParams.SEG_CHANNEL     = SEG_CHANNEL;
handles.nd.segParams.SEG_THRESH      = SEG_THRESH;
handles.nd.segParams.SEG_SMOOTHING   = SEG_SMOOTHING;
handles.nd.segParams.SEG_GAMMA       = SEG_GAMMA;

handles.nd.segParams.useGamma       = get(handles.checkbox_useSegGamma, 'Value');
handles.nd.segParams.useOtsuThresh  = get(handles.checkbox_useOtsuThresh, 'Value');

show_labels = get(handles.checkbox_displaySurfLabels, 'value');

% handles.nd =  segAxis3SIM2(handles.nd, SEG_CHANNEL, SEG_THRESH, SEG_SMOOTHING, []);

if handles.gonad.isConfocal
    
    handles.nd =  segAxis3Confocalv2(handles.nd,...
        SEG_CHANNEL,...
        SEG_THRESH,...
        SEG_SMOOTHING,...
        handles.gonad.MIN_SURF_AREA,...
        handles.gonad.MIN_REG_VOL);
    
else
    
    handles.nd = segAxis3DV(handles.nd);
    
end


handles.nd.surfOn = handles.nd.surfL;
handles.nd.surfOff = 0*handles.nd.surfOn;

% vizSegInitial(handles.nd, handles.nd.surfOn, handles.nd.surfOff, 112, 1, show_labels);
% end


guidata(hObject, handles);

'==============================  FINISHED =============================='


end


% ===========================================================----
%
%  FOCUS-ASSISTED OPTIMIZATION BUTTON
%
% ============================================================--

% --- Executes on button press in button_focusOpti.
function button_focusOpti_Callback(hObject, eventdata, handles)

% ============================================================--
%
%  for now, we assume there are only foci in at most two channels
%
% ============================================================--

'==============================  FOCUS-ASSISTED OPTIMIZING =============================='

nd = handles.nd;

focusChannel = [...
    str2double(get(handles.focus_channel_1, 'String')),...
    str2double(get(handles.focus_channel_2, 'String'))];


nd.focusChannel = focusChannel;

focusAxisTable = get(handles.focusAxisTable, 'data');
nd.focusAxisTable = focusAxisTable;

% make focusData
nd = makeFocusData(nd);

% optimize the segmentation
nd = optiAxis3SmartFast2v3(nd, 0);

% make the final axis mask
nd = makeAxisMask(nd, handles.labelAxesInVolumeOrder);

% display the seg with foci
nd = displayFinalSeg(nd);
handles.nd = nd;
handles = load_axis_mask_status(handles);


% add all possible surfs to the list box
handles.inds = 1:length(handles.nd.ppS);
vals = [];

for ii = 1:length(handles.inds)
    
    inds_str{ii} = num2str(handles.inds(ii));
    
    if isfield(handles.nd, 'optiData')
        if handles.nd.optiData.status(handles.inds(ii))
            vals = [vals, ii];
        end
    else
        if handles.nd.surf_status(handles.inds(ii))
            vals = [vals, ii];
        end
    end
end

set(handles.listbox_surfs, 'String', inds_str, 'Value', vals);

handles.vals = vals;
handles.inds_str = inds_str;


guidata(hObject, handles);


end


% ============================================================--
%
%  ITERATIVE TRACING FROM FOCI BUTTON
%
% ============================================================--


% --- Executes on button press in button_traceFromFoci.
function button_traceFromFoci_Callback(hObject, eventdata, handles)


nd = handles.nd;

stepsPerPixel = 1;
numSteps = 100;

crop = 5;
initTanVec = [1,0,0];

try
    eval(['spots = handles.nd.spots_im' num2str(nd.focusFittingParams.channel) ';']);
catch
    'no spots in selected channel!'
end

eval(['im = handles.nd.im' num2str(handles.gonad.SEG_CHANNEL), ';']);

im = double(im);
mk = ones(size(im));

for ii = 1:length(spots)
    
    initPosition = spots(ii).r;
    
    trace = simpleIterativeTraceAxis3(im, mk, initPosition, initTanVec, crop, stepsPerPixel, numSteps);
    
    spots(ii).iterativeTrace = trace;
    
end

figNum = 90;

display_nucleus_vol_viz(handles, figNum, 0);

ff = figure(figNum);
set(ff, 'Alphamap', (0:64)*.7/64);

vizSpotsConfocal(nd, figNum, 0, [1,0,0], nd.focusFittingParams.channel);

hold on

for ii = 1:length(spots)
    
    plot3(...
        spots(ii).iterativeTrace(:,2),...
        spots(ii).iterativeTrace(:,1),...
        spots(ii).iterativeTrace(:,3),...
        '-', 'LineWidth', 1,...
        'color', [0,1,1]);
    
end

'';


end


% ============================================================--
%
%  FOCUS-FINDING BUTTON
%
% ============================================================--


% --- Executes on button press in button_findFoci.
function button_findFoci_Callback(hObject, eventdata, handles)

nd = handles.nd;

if handles.focus_channel_1_status
    spot_chan =  str2double(get(handles.focus_channel_1, 'String'));
    thresh =     str2double(get(handles.focus_thresh_1, 'String'));
    smoothing =  str2double(get(handles.focus_smoothing_1, 'String'));
elseif handles.focus_channel_2_status
    spot_chan =  str2double(get(handles.focus_channel_2, 'String'));
    thresh =     str2double(get(handles.focus_thresh_2, 'String'));
    smoothing =  str2double(get(handles.focus_smoothing_2, 'String'));
end

try
    eval(['spots = handles.nd.spots_im' num2str(spot_chan) ';']);
catch
    'no spots in selected channel!'
end

eval(['im = handles.nd.im' num2str(spot_chan), ';']);

im = autogain(im);
im = im - mean(im(:));



% find foci

'==================== FINDING FOCI ======================='

REDO = 1;
nd = spotAxis3(nd, spot_chan, thresh, REDO, smoothing);

'==================== FINISHED FINDING FOCI ======================='


% match foci to axes
% nd = spotAxisID3(nd, 10, 3/4);

% vizAxisTraceConfocal(nd, 88, 1, [1,0,0]);
% vizSpotsConfocal(nd, 88, 0, [0,1,0]);

figNum = 89;

try
    cc = handles.nd.colormap;
    vizAxisTraceConfocal(nd, figNum, 1, cc);
    vizSpotsConfocal(nd, figNum, 0, cc, spot_chan);
catch
    cc = [1,1,1];
    vizIm(im, figNum);
    colormap(gray);
    ff = figure(figNum);
    set(ff, 'Alphamap', (0:64)*.3/64);
    vizSpotsConfocal(nd, figNum, 0, cc, spot_chan);
end

handles.nd = nd;

guidata(hObject, handles);

end

% ============================================================--
%
%  SHOW FOCI BUTTON
%
% ============================================================--

% --- Executes on button press in button_showFoci.
function button_showFoci_Callback(hObject, eventdata, handles)

nd = handles.nd;

if handles.focus_channel_1_status
    spot_chan_ind = 1;
    spot_chan     = str2double(get(handles.focus_channel_1,   'String'));
    thresh        = str2double(get(handles.focus_thresh_1,    'String'));
    chan          = str2double(get(handles.focus_channel_1,   'String'));
    smoothing     = str2double(get(handles.focus_smoothing_1, 'String'));
    
elseif handles.focus_channel_2_status
    spot_chan_ind = 2;
    spot_chan     = str2double(get(handles.focus_channel_2,   'String'));
    thresh        = str2double(get(handles.focus_thresh_2,    'String'));
    chan          = str2double(get(handles.focus_channel_2,   'String'));
    smoothing     = str2double(get(handles.focus_smoothing_2, 'String'));
end

try
    eval(['spots = handles.nd.spots_im' num2str(spot_chan) ';']);
catch
    'no spots in selected channel!'
end


figNum = 89;
cc = [1,0,0];

% eval(['im = handles.nd.im' num2str(spot_chan), ';']);
% im = autogain(im);
% im = im - mean(im(:));

% vizIm(im, figNum);
% colormap(gray);

display_nucleus_vol_viz(handles, figNum, 0);

ff = figure(figNum);
set(ff, 'Alphamap', (0:64)*.3/64);

vizSpotsConfocal(nd, figNum, 0, cc, chan);



end


% ============================================================--
%
%  SWITCH FOCUS BUTTON
%
% ============================================================--
% --- Executes on button press in button_switchFoci.
function button_switchFoci_Callback(hObject, eventdata, handles)


if handles.focus_channel_1_status
    spot_chan = str2double(get(handles.focus_channel_1, 'String'));
    spot_chan_ind = 1;
elseif handles.focus_channel_2_status
    spot_chan = str2double(get(handles.focus_channel_2, 'String'));
    spot_chan_ind = 2;
end


focusAxisTable = get(handles.focusAxisTable, 'data');


focusAxisTable([1,2], spot_chan_ind) = focusAxisTable([2,1], spot_chan_ind);

set(handles.focusAxisTable, 'data', focusAxisTable);

handles.nd.focusAxisTable = focusAxisTable;

guidata(hObject, handles);



end

% ============================================================--
%
%  APPLY CHANGES TO FOCUS-AXIS TABLE
%
% ============================================================--


% --- Executes on button press in button_applyChanges.
function button_applyChanges_Callback(hObject, eventdata, handles)


fd = handles.nd.fd;

focusAxisTable = handles.nd.focusAxisTable;

focusAxisTable_new = get(handles.focusAxisTable, 'data');


for nn = 1:2
    chrID = focusAxisTable(:, nn);
    chrID_new = focusAxisTable_new(:, nn);
    chan = handles.nd.focusChannel(nn);
    
    for ii = 1:length(chrID)
        for jj = 1:length(fd)
            
            if fd{jj}.chrID==chrID(ii) && fd{jj}.channel==chan
                
                fd{jj}.chrID = chrID_new(ii);
                
            end
        end
    end
end


handles.nd.focusAxisTable = focusAxisTable_new;
handles.nd.fd =             fd;
handles.nd =                displayFinalSeg(handles.nd);

guidata(hObject, handles);

end

% ============================================================--
%
%  FIND MIN PATH BETWEEN FOCI
%
% ============================================================--

% --- Executes on button press in button_findMinPath.
function button_findMinPath_Callback(hObject, eventdata, handles)


scale = [str2double(get(handles.im1_scale,'String')),...
    str2double(get(handles.im2_scale,'String')),...
    str2double(get(handles.im3_scale,'String')),...
    str2double(get(handles.im4_scale,'String'))];

offset = [str2double(get(handles.im1_off,'String')),...
    str2double(get(handles.im2_off,'String')),...
    str2double(get(handles.im3_off,'String')),...
    str2double(get(handles.im4_off,'String'))];

gamma = [str2double(get(handles.im1_gamma,'String')),...
    str2double(get(handles.im2_gamma,'String')),...
    str2double(get(handles.im3_gamma,'String')),...
    str2double(get(handles.im4_gamma,'String'))];

handles.nd.displaySettings.scale = scale;
handles.nd.displaySettings.offset = offset;
handles.nd.displaySettings.gamma = gamma;


nd = handles.nd;

focusChannel = [...
    str2double(get(handles.focus_channel_1, 'String')),...
    str2double(get(handles.focus_channel_2, 'String'))];

nd.focusChannel = focusChannel;

focusAxisTable = get(handles.focusAxisTable, 'data');
nd.focusAxisTable = focusAxisTable;

% make focusData
nd = makeFocusData(nd);

% calc min path
nd = calcMinPathBetweenFoci(nd);


figNum = 899;
display_nucleus_vol_viz(handles, figNum, 0);

ff = figure(figNum);
set(ff, 'Alphamap', (0:64)*.3/64);

vizVolume3(nd.pathData{1}.pathImage, [1,0,0], figNum, 0, 0, 1);


'';

end




% ============================================================--
%
%  OPTIMIZATION BUTTONS
%
% ============================================================--

% --- Executes on button press in button_optimize.
function button_optimize_Callback(hObject, eventdata, handles)

handles.resave_nd_flag = 1;

handles.nd = optiAxis3SmartFast2v3(handles.nd);
handles.nd = makeAxisMask(handles.nd, handles.labelAxesInVolumeOrder);
handles.nd = displayFinalSeg(handles.nd);
handles = load_axis_mask_status(handles);

% add all possible surfs to the list box
handles.inds = 1:length(handles.nd.ppS);
vals = [];

for ii = 1:length(handles.inds)
    
    inds_str{ii} = num2str(handles.inds(ii));
    
    if isfield(handles.nd, 'optiData')
        
        if handles.nd.optiData.status(handles.inds(ii))
            vals = [vals, ii];
        end
        
    else
        
        if handles.nd.surf_status(handles.inds(ii))
            vals = [vals, ii];
        end
        
    end
    
    
end

set(handles.listbox_surfs, 'String', inds_str, 'Value', vals);

handles.vals = vals;
handles.inds_str = inds_str;

guidata(hObject, handles);

end


% --- Executes on button press in button_optimize_old.
function button_optimize_old_Callback(hObject, eventdata, handles)

handles.resave_nd_flag = 1;

if (get(handles.radio_intensity,'Value') == get(handles.radio_intensity,'Max'))
    order_param_str = 'intensity';
end
if (get(handles.radio_orientation,'Value') == get(handles.radio_orientation,'Max'))
    order_param_str = 'orientation';
end
if (get(handles.radio_random,'Value') == get(handles.radio_random,'Max'))
    order_param_str = 'random';
end

% MAX_AXIS_LENGTH = str2double(get(handles.seg_gamma, 'String'));

MAX_AXIS_LENGTH = 999;

handles.nd = optiAxis3Quick(handles.nd, order_param_str, MAX_AXIS_LENGTH);

handles.nd = displayFinalSeg(handles.nd);

% add all possible surfs to the list box
handles.inds = 1:length(handles.nd.ppS);
vals = [];

for ii = 1:length(handles.inds)
    
    inds_str{ii} = num2str(handles.inds(ii));
    
    if handles.nd.surf_status(handles.inds(ii))
        vals = [vals, ii];
    end
    
end

set(handles.listbox_surfs, 'String', inds_str, 'Value', vals);

handles.vals = vals;
handles.inds_str = inds_str;

guidata(hObject, handles);

end

% ============================================================--
%
%  *** NUCLEUS STATUS BUTTONS ***
%
% ============================================================--


% --- Executes on button press in button_ignore.
function button_finalize_Callback(hObject, eventdata, handles)

ID_NUMBER = 1;

if handles.gonad.finalized_nukes(handles.index_selected)==ID_NUMBER
    handles.gonad.finalized_nukes(handles.index_selected) = 0;
    set(handles.button_finalize, 'FontWeight', 'Normal');
else
    handles.gonad.finalized_nukes(handles.index_selected) = ID_NUMBER;
    set(handles.button_finalize, 'FontWeight', 'Bold');
    set(handles.button_ignore, 'FontWeight', 'Normal');
    set(handles.button_error, 'FontWeight', 'Normal');
end

if handles.useGonadMaskProjection
    handles = update_nuke_mask(handles, ID_NUMBER);
end

display_gonad(handles, hObject);

gonad = handles.gonad;
save( [handles.dd filesep 'gonad.mat'], 'gonad');

if handles.resave_nd_flag
    
    nd = handles.nd;
    
    dd = [handles.dd filesep ...
        handles.nukeDir filesep ...
        num2str(handles.gonad.good_nukes(handles.index_selected)) filesep...
        'nd.mat']
    
    save(dd, 'nd');
    
    'SAVED'
    
end

set(handles.textbox_number_finalize, 'string', num2str(sum(handles.gonad.finalized_nukes==1)));
set(handles.textbox_number_error, 'string', num2str(sum(handles.gonad.finalized_nukes==2)));
set(handles.textbox_number_ignore, 'string', num2str(sum(handles.gonad.finalized_nukes==-1)));

guidata(hObject, handles);

end

% --- Executes on button press in button_finalize.
function button_ignore_Callback(hObject, eventdata, handles)

ID_NUMBER = -1;

if handles.gonad.finalized_nukes(handles.index_selected)==ID_NUMBER
    handles.gonad.finalized_nukes(handles.index_selected) = 0;
    set(handles.button_ignore, 'FontWeight', 'Normal');
else
    handles.gonad.finalized_nukes(handles.index_selected) = ID_NUMBER;
    set(handles.button_ignore, 'FontWeight', 'Bold');
    set(handles.button_finalize, 'FontWeight', 'Normal');
    set(handles.button_error, 'FontWeight', 'Normal');
end

if handles.useGonadMaskProjection
    handles = update_nuke_mask(handles, ID_NUMBER);
end

display_gonad(handles, hObject);

gonad = handles.gonad;
save( [handles.dd filesep 'gonad.mat'], 'gonad');

set(handles.textbox_number_finalize, 'string', num2str(sum(handles.gonad.finalized_nukes==1)));
set(handles.textbox_number_error, 'string', num2str(sum(handles.gonad.finalized_nukes==2)));
set(handles.textbox_number_ignore, 'string', num2str(sum(handles.gonad.finalized_nukes==-1)));

guidata(hObject, handles);

end

% --- Executes on button press in button_error.
function button_error_Callback(hObject, eventdata, handles)

ID_NUMBER = 2;

if handles.gonad.finalized_nukes(handles.index_selected)==ID_NUMBER
    handles.gonad.finalized_nukes(handles.index_selected) = 0;
    set(handles.button_error, 'FontWeight', 'Normal');
else
    handles.gonad.finalized_nukes(handles.index_selected) = ID_NUMBER;
    set(handles.button_error, 'FontWeight', 'Bold');
    set(handles.button_ignore, 'FontWeight', 'Normal');
    set(handles.button_finalize, 'FontWeight', 'Normal');
end

if handles.useGonadMaskProjection
    handles = update_nuke_mask(handles, ID_NUMBER);
end

display_gonad(handles, hObject);

gonad = handles.gonad;
save( [handles.dd filesep 'gonad.mat'], 'gonad');

set(handles.textbox_number_finalize, 'string', num2str(sum(handles.gonad.finalized_nukes==1)));
set(handles.textbox_number_error, 'string', num2str(sum(handles.gonad.finalized_nukes==2)));
set(handles.textbox_number_ignore, 'string', num2str(sum(handles.gonad.finalized_nukes==-1)));

guidata(hObject, handles);

end

% ============================================================--
%
%   GONAD MASK UPDATE FUNCTION
%
% ============================================================--


function handles = update_nuke_mask(handles, ID_NUMBER)

% good nucleus
if ID_NUMBER==1
    cc = [55,255,55];
    
    % ignore
elseif ID_NUMBER==-1
    cc = [155,155,155];
    
    % segmentation error
elseif ID_NUMBER==2
    cc = [255,55,0];
    
    % should never get here
else
    'should not get to this line!'
end

if handles.gonad.finalized_nukes(handles.index_selected)==ID_NUMBER
    
    imb = handles.mkFC_projection;
    mk = max(handles.mkF==handles.gonad.good_nukes(handles.index_selected),[],3);
    
    imR = imb(:,:,1);
    imG = imb(:,:,2);
    imB = imb(:,:,3);
    
    imR(mk) = cc(1);
    imG(mk) = cc(2);
    imB(mk) = cc(3);
    
    handles.mkFC_projection(:,:,1) = imR;
    handles.mkFC_projection(:,:,2) = imG;
    handles.mkFC_projection(:,:,3) = imB;
    
else
    
    mk = max(handles.mkF==handles.gonad.good_nukes(handles.index_selected),[],3);
    
    imb = handles.mkFC_projection;
    imR = imb(:,:,1);
    imG = imb(:,:,2);
    imB = imb(:,:,3);
    
    imb_orig = handles.mkFC_projection_orig;
    
    imR_orig = imb_orig(:,:,1);
    imG_orig = imb_orig(:,:,2);
    imB_orig = imb_orig(:,:,3);
    
    imR(mk) = imR_orig(mk);
    imG(mk) = imG_orig(mk);
    imB(mk) = imB_orig(mk);
    
    handles.mkFC_projection(:,:,1) = imR;
    handles.mkFC_projection(:,:,2) = imG;
    handles.mkFC_projection(:,:,3) = imB;
    
end


end

function handles = load_axis_mask_status(handles)

numAxes = size(handles.nd.colormap,1);

for jj = 1:length(handles.nd.groups_final)
    axisAreas(jj) = sum([handles.nd.ppR(handles.nd.groups_final(jj).regs).Area]);
end

[ss, ord] = sort(axisAreas, 'Descend');

for ii = 1:10
    
    if ii > numAxes
        eval(['set(handles.button_includeAxis_' num2str(ii) ', ''BackgroundColor'', [1,1,1]/3);']);
    else
        eval(['set(handles.button_includeAxis_' num2str(ii) ', ''BackgroundColor'', handles.nd.colormap(ii,:));']);
    end
    
end

if ~isfield(handles.nd, 'mkL_final_includeFlags')
    handles.nd.mkL_final_includeFlags = 0*(1:numAxes);
end

flags = handles.nd.mkL_final_includeFlags;

if length(flags) < numAxes
    flags(end+1:end+(numAxes - length(flags))) = 0;
end

handles.nd.mkL_final_includeFlags = flags;

for ii = 1:10
    
    if ii > length(flags)
        eval(['set(handles.button_includeAxis_' num2str(ii) ', ''String'', ''Blank'');']);
        eval(['set(handles.button_includeAxis_' num2str(ii) ', ''FontWeight'', ''Normal'');']);
    else
        if flags(ii)
            eval(['set(handles.button_includeAxis_' num2str(ii) ', ''String'', ''Yes'');']);
            eval(['set(handles.button_includeAxis_' num2str(ii) ', ''FontWeight'', ''Bold'');']);
        else
            eval(['set(handles.button_includeAxis_' num2str(ii) ', ''String'', ''No'');']);
            eval(['set(handles.button_includeAxis_' num2str(ii) ', ''FontWeight'', ''Normal'');']);
        end
    end
    
    
end
end


% ============================================================--
%
%  ***SAVE NUCLEUS BUTTON***
%
% ============================================================--


% --- Executes on button press in button_saveNucleus.
function button_saveNucleus_Callback(hObject, eventdata, handles)

nd = handles.nd;

save(...
    [handles.dd filesep ...
    handles.nukeDir filesep ...
    num2str(handles.gonad.good_nukes(handles.index_selected)) filesep...
    'nd.mat'],...
    'nd');


fname = [handles.dd filesep ...
    handles.nukeDir filesep ...
    num2str(handles.gonad.good_nukes(handles.index_selected)) filesep...
    'nd.mat']

figure(1);

'==============================  SAVED =============================='

end

% ============================================================--
%
%  SEGMENTATION SURFACES -- DISPLAY/EDIT
%
% ============================================================--


% --- Executes on button press in checkbox_displaySurfLabels.
function checkbox_displaySurfLabels_Callback(hObject, eventdata, handles)

% show_labels = get(handles.checkbox_displaySurfLabels, 'value');
%
% figure(11);
% pos = get(gca, 'CameraPosition');
% vizSegInitial(handles.nd, 11, 1, show_labels);
% set(gca, 'CameraPosition', pos);

end

% --- Executes on button press in button_toggleSurfLabels.
function button_toggleSurfLabels_Callback(hObject, eventdata, handles)
end


% --- Executes on selection change in listbox_surfs.
function listbox_surfs_Callback(hObject, eventdata, handles)

handles.resave_nd_flag = 1;

val = get(handles.listbox_surfs, 'Value');

changed_val = setxor(handles.vals, val);

if isempty(intersect(handles.vals, changed_val))
    handles.vals = [handles.vals, changed_val];
else
    handles.vals = setdiff(handles.vals, changed_val);
end

new_ind = handles.inds(changed_val);

if isfield(handles.nd, 'optiData')
    handles.nd.optiData.status(new_ind) = ~handles.nd.optiData.status(new_ind);
else
    handles.nd.surf_status(new_ind) = ~handles.nd.surf_status(new_ind);
end

handles.nd = makeAxisMask(handles.nd, handles.labelAxesInVolumeOrder);

figNum = 111;
handles.nd = displayFinalSeg(handles.nd, figNum);

handles.nd = rmfield(handles.nd,'mkL_final_includeFlags');

handles = load_axis_mask_status(handles);

guidata(hObject, handles);


end


function handles = display_final_seg(handles)

% handles.nd = makeAxisMask(handles.nd, handles.labelAxesInVolumeOrder);
%
% cc = colormap(lines(max(handles.nd.mkL_final(:))));
% handles.cc = cc;
%
% vizSegColor(handles.nd.mkL_final, 111, 1, cc);
% set(gcf, 'position', [420,1040,560,420]);

% vizVolume2Color(handles.nd.skL_final, 0, 1111, 1, 1, cc);
% set(gcf, 'position', [813,1065,560,420]);

% pp = regionprops(handles.nd.skL_final,'Area');
% sort([pp.Area])

end

% ============================================================--
%
%  CONVENTIONAL TRACING BUTTON
%
% ============================================================--

% --- Executes on button press in button_trace.
function button_trace_Callback(hObject, eventdata, handles)


try
    vizAxisTraceConfocal(handles.nd, 1111, 1, handles.nd.colormap);
catch
    
    if ~isfield(handles.nd, 'traceParams')
        handles.nd.traceParams.CROP = 4;
        handles.nd.traceParams.INTERPS_PER_PIXEL = 5;
        handles.nd.traceParams.STEPS_PER_PIXEL = 1;
        handles.nd.traceParams.CONV_LENGTH = 3;
    end
    
    handles.nd = makeAxisDataConfocal(handles.nd);
    handles.nd = traceSkel5(handles.nd);
    guidata(hObject, handles);
    
    vizAxisTraceConfocal(handles.nd, 1111, 1, handles.nd.colormap);
    
end


display_axis_images(handles);



guidata(hObject, handles);


end



% ============================================================--
%
%  SHOW SURFACE SCORES BUTTON
%
% ============================================================--


% --- Executes on button press in button_showScores.
function button_showScores_Callback(hObject, eventdata, handles)

nd = handles.nd;

% load calculated score parameters and probability parameters
[scoreParams, probGoodParams] = loadScoreFunctionParams3;

% load the data set used to calculate those parameters
% i.e., the training data set
dataGlobal = loadGlobalData3;

% make surface properties for this nucleus
[props, status] = makePropsArray2v3(nd);
data(1).props = props;
data(1).status = status;

% calculate surface properties for this nucleus using surface data
% and certain averages from the training data set
avProps = calcAvSurfProps3(data, dataGlobal);

% now we can calculate scores
scores = calcSurfScoreFunction3(avProps, scoreParams);

% and the probabilities that each surface is good
probGood = calcSurfGoodProb(scores, probGoodParams);

vizSegScores(handles.nd, 113, 1, probGood);


end


% ============================================================--
%
%  PER-AXIS 3D IMAGE DISPLAY
%
% ============================================================--


function display_axis_images(handles)

NUM_AXES = length(handles.nd.sdata);

for ii = 1:min(6,NUM_AXES)
    
    [im, channels] = make_axis_images(handles, ii);
    
    figure(100+ii);clf;
    
    if sum(~~channels)==1 && get(handles.checkbox_use_colormap, 'Value')
        
        im = squeeze(im(:,:,:,channels(find(~~channels))));
        
        vol3d('CData', autogain((double(im)+1)),'texture','3D','AlphaData', im);
        
        cc = flipud(colormap(jet(255)));
        colormap(cc);
        
    else
        
        vol3d('CData', im,'texture','3D','AlphaData', im);
        
    end
    
    daspect([1 1 1]);
    view(3);
    axis tight;
    %  camlight left;
    %  camlight right;
    axis vis3d;
    %  lighting none; %without this, no intensity mapping
    lighting flat;
    material dull;
    
    rotate3d;
    set(gca, 'Position', [0,0,1,1]);
    
    ff = gcf;
    set(gcf, 'color', [0,0,0]);
    
    handles.alphamap = str2double(get(handles.alpha, 'String'));
    set(ff, 'Alphamap', (0:64)*handles.alphamap/64);
    
end

end


function [im, channels] = make_axis_images(handles, axis_index)


scale = [str2double(get(handles.im1_scale,'String')),...
    str2double(get(handles.im2_scale,'String')),...
    str2double(get(handles.im3_scale,'String')),...
    str2double(get(handles.im4_scale,'String'))];

offset = [str2double(get(handles.im1_off,'String')),...
    str2double(get(handles.im2_off,'String')),...
    str2double(get(handles.im3_off,'String')),...
    str2double(get(handles.im4_off,'String'))];

gamma = [str2double(get(handles.im1_gamma,'String')),...
    str2double(get(handles.im2_gamma,'String')),...
    str2double(get(handles.im3_gamma,'String')),...
    str2double(get(handles.im4_gamma,'String'))];

rr = ceil(handles.nd.sdata(axis_index).rr);
cc = ceil(handles.nd.sdata(axis_index).cc);
zz = ceil(handles.nd.sdata(axis_index).zz);

try
    im1 = handles.nd.sdata(axis_index).im{1};
    im2 = handles.nd.sdata(axis_index).im{2};
    im3 = handles.nd.sdata(axis_index).im{3};
    im4 = handles.nd.sdata(axis_index).im{4};
catch
    im1 = handles.im1(rr,cc,zz);
    im2 = handles.im2(rr,cc,zz);
    im3 = handles.im3(rr,cc,zz);
    im4 = handles.im4(rr,cc,zz);
end

try
    mk = uint8(handles.nd.sdata(axis_index).mk);
catch
    mk = uint8(handles.nd.mkL_final(rr,cc,zz)==handles.nd.sdata(axis_index).label);
end

im1 = (scale(1)*(im1 - mean(im1(:))*offset(1)));
im2 = (scale(2)*(im2 - mean(im2(:))*offset(2)));
im3 = (scale(3)*(im3 - mean(im3(:))*offset(3)));
im4 = (scale(4)*(im4 - mean(im4(:))*offset(4)));

im1  = autogain(double(im1).^(gamma(1)));
im2  = autogain(double(im2).^(gamma(2)));
im3  = autogain(double(im3).^(gamma(3)));
im4  = autogain(double(im4).^(gamma(4)));

ord = [handles.im1_index, handles.im2_index, handles.im3_index, handles.im4_index];

channels = (ord);

r_inds = find(ord==1);
g_inds = find(ord==2);
b_inds = find(ord==3);

imr = autogain(0*im1);
img = autogain(0*im1);
imb = autogain(0*im1);

% for ii = r_inds
%     imr = imr + eval(['im' num2str(ii) ';']);
% end
% for ii = g_inds
%     img = img + eval(['im' num2str(ii) ';']);
% end
% for ii = b_inds
%     imb = imb + eval(['im' num2str(ii) ';']);
% end

if sum(handles.im1_index==1)
    imr = imr + im1;
end
if sum(handles.im1_index==2)
    img = img + im1;
end
if sum(handles.im1_index==3)
    imb = imb + im1;
end
if sum(handles.im2_index==1)
    imr = imr + im2;
end
if sum(handles.im2_index==2)
    img = img + im2;
end
if sum(handles.im2_index==3)
    imb = imb + im2;
end
if sum(handles.im3_index==1)
    imr = imr + im3;
end
if sum(handles.im3_index==2)
    img = img + im3;
end
if sum(handles.im3_index==3)
    imb = imb + im3;
end
if sum(handles.im4_index==1)
    imr = imr + im4;
end
if sum(handles.im4_index==2)
    img = img + im4;
end
if sum(handles.im4_index==3)
    imb = imb + im4;
end

if get(handles.checkbox_use_smoothing,'Value')
    imr = imresize3(imr, 2);
    img = imresize3(img, 2);
    imb = imresize3(imb, 2);
    mk = imresize3(mk, 2);
end

if get(handles.checkbox_useMask, 'Value')
    imr = imr.*mk;
    img = img.*mk;
    imb = imb.*mk;
end

rr = [255, 0, 60]/255;
gg = [60, 255, 0]/255;
bb = [0, 125, 255]/255;

rr = [200, 25, 25]/255;
gg = [25, 200, 25]/255;
bb = [25, 25, 200]/255;

% rr = [1, 1, 0];
% gg = [1, 0, 1];
% bb = [0, 1, 1];
%
% rr = [1, 1/2, 1/2];
% gg = [1/2, 1, 1/2];
% bb = [1/2, 1/2, 1];

imR = imr*rr(1) + img*gg(1) + imb*bb(1);
imG =  imr*rr(2) + img*gg(2) + imb*bb(2);
imB = imr*rr(3) + img*gg(3) + imb*bb(3);

im = cat(4,...
    imr*rr(1) + img*gg(1) + imb*bb(1),...
    imr*rr(2) + img*gg(2) + imb*bb(2),...
    imr*rr(3) + img*gg(3) + imb*bb(3));

end



function alpha_Callback(hObject, eventdata, handles)

handles.alphamap = str2double(get(handles.alpha, 'String'));
set(gcf, 'Alphamap', (0:64)*handles.alphamap/64);
guidata(hObject, handles);

end


% --- Executes on button press in checkbox_use_colormap.
function checkbox_use_colormap_Callback(hObject, eventdata, handles)

display_nucleus_vol_viz(handles);

end


% ============================================================--
%
%  FOCUS CHANNEL SELECTION BUTTONS
%
% ============================================================--


% --- Executes on button press in button_focus_channel_1.
function button_focus_channel_1_Callback(hObject, eventdata, handles)

handles.focus_channel_1_status = ~handles.focus_channel_1_status;
handles.focus_channel_2_status = ~handles.focus_channel_2_status;

if handles.focus_channel_1_status
    set(handles.button_focus_channel_1, 'FontWeight', 'Bold');
else
    set(handles.button_focus_channel_1, 'FontWeight', 'Normal');
end
if handles.focus_channel_2_status
    set(handles.button_focus_channel_2, 'FontWeight', 'Bold');
else
    set(handles.button_focus_channel_2, 'FontWeight', 'Normal');
end

guidata(hObject, handles);
end


% --- Executes on button press in button_focus_channel_2.
function button_focus_channel_2_Callback(hObject, eventdata, handles)

handles.focus_channel_1_status = ~handles.focus_channel_1_status;
handles.focus_channel_2_status = ~handles.focus_channel_2_status;


if handles.focus_channel_1_status
    set(handles.button_focus_channel_1, 'FontWeight', 'Bold');
else
    set(handles.button_focus_channel_1, 'FontWeight', 'Normal');
end
if handles.focus_channel_2_status
    set(handles.button_focus_channel_2, 'FontWeight', 'Bold');
else
    set(handles.button_focus_channel_2, 'FontWeight', 'Normal');
end
guidata(hObject, handles);

end


% ============================================================--
%
%  AXIS INCLUSION BUTTONS
%
% ============================================================--

% --- Executes on button press in button_includeAxis_all.
function button_includeAxis_all_Callback(hObject, eventdata, handles)


handles.resave_nd_flag = 1;

inds = 1:6;

for ind = inds

    if handles.nd.mkL_final_includeFlags(ind)
        handles.nd.mkL_final_includeFlags(ind) = 0;
        eval(['set(handles.button_includeAxis_' num2str(ind) ', ''String'', ''No'');']);
        eval(['set(handles.button_includeAxis_' num2str(ind) ', ''FontWeight'', ''Normal'');']);
    else
        handles.nd.mkL_final_includeFlags(ind) = 1;
        eval(['set(handles.button_includeAxis_' num2str(ind) ', ''String'', ''Yes'');']);
        eval(['set(handles.button_includeAxis_' num2str(ind) ', ''FontWeight'', ''Bold'');']);
    end
end

guidata(hObject, handles);


end


% --- Executes on button press in button_includeAxis_1.
function button_includeAxis_1_Callback(hObject, eventdata, handles)

handles.resave_nd_flag = 1;

ind = 1;

if handles.nd.mkL_final_includeFlags(ind)
    handles.nd.mkL_final_includeFlags(ind) = 0;
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''String'', ''No'');']);
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''FontWeight'', ''Normal'');']);
else
    handles.nd.mkL_final_includeFlags(ind) = 1;
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''String'', ''Yes'');']);
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''FontWeight'', ''Bold'');']);
end
guidata(hObject, handles);
end

% --- Executes on button press in button_includeAxis_2.
function button_includeAxis_2_Callback(hObject, eventdata, handles)

handles.resave_nd_flag = 1;

ind = 2;

if handles.nd.mkL_final_includeFlags(ind)
    handles.nd.mkL_final_includeFlags(ind) = 0;
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''String'', ''No'');']);
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''FontWeight'', ''Normal'');']);
else
    handles.nd.mkL_final_includeFlags(ind) = 1;
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''String'', ''Yes'');']);
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''FontWeight'', ''Bold'');']);
end
guidata(hObject, handles);
end

% --- Executes on button press in button_includeAxis_3.
function button_includeAxis_3_Callback(hObject, eventdata, handles)
handles.resave_nd_flag = 1;

ind = 3;

if handles.nd.mkL_final_includeFlags(ind)
    handles.nd.mkL_final_includeFlags(ind) = 0;
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''String'', ''No'');']);
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''FontWeight'', ''Normal'');']);
else
    handles.nd.mkL_final_includeFlags(ind) = 1;
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''String'', ''Yes'');']);
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''FontWeight'', ''Bold'');']);
end
guidata(hObject, handles);
end

% --- Executes on button press in button_includeAxis_4.
function button_includeAxis_4_Callback(hObject, eventdata, handles)

handles.resave_nd_flag = 1;

ind = 4;

if handles.nd.mkL_final_includeFlags(ind)
    handles.nd.mkL_final_includeFlags(ind) = 0;
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''String'', ''No'');']);
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''FontWeight'', ''Normal'');']);
else
    handles.nd.mkL_final_includeFlags(ind) = 1;
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''String'', ''Yes'');']);
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''FontWeight'', ''Bold'');']);
end
guidata(hObject, handles);
end

% --- Executes on button press in button_includeAxis_5.
function button_includeAxis_5_Callback(hObject, eventdata, handles)

handles.resave_nd_flag = 1;

ind = 5;

if handles.nd.mkL_final_includeFlags(ind)
    handles.nd.mkL_final_includeFlags(ind) = 0;
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''String'', ''No'');']);
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''FontWeight'', ''Normal'');']);
else
    handles.nd.mkL_final_includeFlags(ind) = 1;
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''String'', ''Yes'');']);
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''FontWeight'', ''Bold'');']);
end
guidata(hObject, handles);
end

% --- Executes on button press in button_includeAxis_6.
function button_includeAxis_6_Callback(hObject, eventdata, handles)

handles.resave_nd_flag = 1;

ind = 6;

if handles.nd.mkL_final_includeFlags(ind)
    handles.nd.mkL_final_includeFlags(ind) = 0;
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''String'', ''No'');']);
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''FontWeight'', ''Normal'');']);
else
    handles.nd.mkL_final_includeFlags(ind) = 1;
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''String'', ''Yes'');']);
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''FontWeight'', ''Bold'');']);
end
guidata(hObject, handles);
end

% --- Executes on button press in button_includeAxis_7.
function button_includeAxis_7_Callback(hObject, eventdata, handles)

handles.resave_nd_flag = 1;

ind = 7;

if handles.nd.mkL_final_includeFlags(ind)
    handles.nd.mkL_final_includeFlags(ind) = 0;
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''String'', ''No'');']);
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''FontWeight'', ''Normal'');']);
else
    handles.nd.mkL_final_includeFlags(ind) = 1;
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''String'', ''Yes'');']);
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''FontWeight'', ''Bold'');']);
end
guidata(hObject, handles);
end

% --- Executes on button press in button_includeAxis_8.
function button_includeAxis_8_Callback(hObject, eventdata, handles)

handles.resave_nd_flag = 1;

ind = 8;

if handles.nd.mkL_final_includeFlags(ind)
    handles.nd.mkL_final_includeFlags(ind) = 0;
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''String'', ''No'');']);
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''FontWeight'', ''Normal'');']);
else
    handles.nd.mkL_final_includeFlags(ind) = 1;
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''String'', ''Yes'');']);
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''FontWeight'', ''Bold'');']);
end
guidata(hObject, handles);
end

% --- Executes on button press in button_includeAxis_9.
function button_includeAxis_9_Callback(hObject, eventdata, handles)
handles.resave_nd_flag = 1;

ind = 9;

if handles.nd.mkL_final_includeFlags(ind)
    handles.nd.mkL_final_includeFlags(ind) = 0;
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''String'', ''No'');']);
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''FontWeight'', ''Normal'');']);
else
    handles.nd.mkL_final_includeFlags(ind) = 1;
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''String'', ''Yes'');']);
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''FontWeight'', ''Bold'');']);
end
guidata(hObject, handles);
end

% --- Executes on button press in button_includeAxis_10.
function button_includeAxis_10_Callback(hObject, eventdata, handles)
handles.resave_nd_flag = 1;

ind = 10;

if handles.nd.mkL_final_includeFlags(ind)
    handles.nd.mkL_final_includeFlags(ind) = 0;
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''String'', ''No'');']);
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''FontWeight'', ''Normal'');']);
else
    handles.nd.mkL_final_includeFlags(ind) = 1;
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''String'', ''Yes'');']);
    eval(['set(handles.button_includeAxis_' num2str(ind) ', ''FontWeight'', ''Bold'');']);
end
guidata(hObject, handles);
end





% ============================================================--
%
%  COLOR CHANNEL BUTTONS
%
% ============================================================--

% --- Executes on button press in button_im1_B.
function button_im1_B_Callback(hObject, eventdata, handles)
% hObject    handle to button_im1_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if sum(handles.im1_index==3)
    handles.im1_index = setxor(handles.im1_index, 3);
    %     set(handles.button_im1_R, 'FontWeight', 'Normal');
    %     set(handles.button_im1_G, 'FontWeight', 'Normal');
    set(handles.button_im1_B, 'FontWeight', 'Normal');
else
    handles.im1_index = [handles.im1_index, 3];
    %     set(handles.button_im1_R, 'FontWeight', 'Normal');
    %     set(handles.button_im1_G, 'FontWeight', 'Normal');
    set(handles.button_im1_B, 'FontWeight', 'Bold');
end

guidata(hObject, handles);
display_gonad(handles, hObject);
end

% --- Executes on button press in button_im1_G.
function button_im1_G_Callback(hObject, eventdata, handles)
% hObject    handle to button_im1_G (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if sum(handles.im1_index==2)
    handles.im1_index = setxor(handles.im1_index, 2);
    %     set(handles.button_im1_R, 'FontWeight', 'Normal');
    set(handles.button_im1_G, 'FontWeight', 'Normal');
    %     set(handles.button_im1_B, 'FontWeight', 'Normal');
else
    handles.im1_index = [handles.im1_index, 2];
    %     set(handles.button_im1_R, 'FontWeight', 'Normal');
    set(handles.button_im1_G, 'FontWeight', 'Bold');
    %     set(handles.button_im1_B, 'FontWeight', 'Normal');
end

guidata(hObject, handles);
display_gonad(handles, hObject);
end

% --- Executes on button press in button_im1_R.
function button_im1_R_Callback(hObject, eventdata, handles)
% hObject    handle to button_im1_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if sum(handles.im1_index==1)
    handles.im1_index = setxor(handles.im1_index, 1);
    set(handles.button_im1_R, 'FontWeight', 'Normal');
    %     set(handles.button_im1_G, 'FontWeight', 'Normal');
    %     set(handles.button_im1_B, 'FontWeight', 'Normal');
else
    handles.im1_index = [handles.im1_index, 1];
    set(handles.button_im1_R, 'FontWeight', 'Bold');
    %     set(handles.button_im1_G, 'FontWeight', 'Normal');
    %     set(handles.button_im1_B, 'FontWeight', 'Normal');
end

guidata(hObject, handles);
display_gonad(handles, hObject);
end

% --- Executes on button press in button_im2_B.
function button_im2_B_Callback(hObject, eventdata, handles)
% hObject    handle to button_im2_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if sum(handles.im2_index==3)
    handles.im2_index = setxor(handles.im2_index, 3);
    %     set(handles.button_im2_R, 'FontWeight', 'Normal');
    %     set(handles.button_im2_G, 'FontWeight', 'Normal');
    set(handles.button_im2_B, 'FontWeight', 'Normal');
else
    handles.im2_index = [handles.im2_index, 3];
    %     set(handles.button_im2_R, 'FontWeight', 'Normal');
    %     set(handles.button_im2_G, 'FontWeight', 'Normal');
    set(handles.button_im2_B, 'FontWeight', 'Bold');
end

guidata(hObject, handles);
display_gonad(handles, hObject);
end

% --- Executes on button press in button_im2_G.
function button_im2_G_Callback(hObject, eventdata, handles)
% hObject    handle to button_im2_G (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if sum(handles.im2_index==2)
    handles.im2_index = setxor(handles.im2_index, 2);
    %     set(handles.button_im2_R, 'FontWeight', 'Normal');
    set(handles.button_im2_G, 'FontWeight', 'Normal');
    %     set(handles.button_im2_B, 'FontWeight', 'Normal');
else
    handles.im2_index = [handles.im2_index, 2];
    %     set(handles.button_im2_R, 'FontWeight', 'Normal');
    set(handles.button_im2_G, 'FontWeight', 'Bold');
    %     set(handles.button_im2_B, 'FontWeight', 'Normal');
end

guidata(hObject, handles);
display_gonad(handles, hObject);
end

% --- Executes on button press in button_im2_R.
function button_im2_R_Callback(hObject, eventdata, handles)
% hObject    handle to button_im2_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if sum(handles.im2_index==1)
    handles.im2_index = setxor(handles.im2_index, 1);
    set(handles.button_im2_R, 'FontWeight', 'Normal');
    %     set(handles.button_im2_G, 'FontWeight', 'Normal');
    %     set(handles.button_im2_B, 'FontWeight', 'Normal');
else
    handles.im2_index = [handles.im2_index, 1];
    set(handles.button_im2_R, 'FontWeight', 'Bold');
    %     set(handles.button_im2_G, 'FontWeight', 'Normal');
    %     set(handles.button_im2_B, 'FontWeight', 'Normal');
end
guidata(hObject, handles);
display_gonad(handles, hObject);
end


% --- Executes on button press in button_im3_B.
function button_im3_B_Callback(hObject, eventdata, handles)
% hObject    handle to button_im3_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if sum(handles.im3_index==3)
    handles.im3_index = setxor(handles.im3_index, 3);
    %     set(handles.button_im3_R, 'FontWeight', 'Normal');
    %     set(handles.button_im3_G, 'FontWeight', 'Normal');
    set(handles.button_im3_B, 'FontWeight', 'Normal');
else
    handles.im3_index = [handles.im3_index, 3];
    %     set(handles.button_im3_R, 'FontWeight', 'Normal');
    %     set(handles.button_im3_G, 'FontWeight', 'Normal');
    set(handles.button_im3_B, 'FontWeight', 'Bold');
end

guidata(hObject, handles);
display_gonad(handles, hObject);
end

% --- Executes on button press in button_im3_G.
function button_im3_G_Callback(hObject, eventdata, handles)
% hObject    handle to button_im3_G (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if sum(handles.im3_index==2)
    handles.im3_index = setxor(handles.im3_index, 2);
    %     set(handles.button_im3_R, 'FontWeight', 'Normal');
    set(handles.button_im3_G, 'FontWeight', 'Normal');
    %     set(handles.button_im3_B, 'FontWeight', 'Normal');
else
    handles.im3_index = [handles.im3_index, 2];
    %     set(handles.button_im3_R, 'FontWeight', 'Normal');
    set(handles.button_im3_G, 'FontWeight', 'Bold');
    %     set(handles.button_im3_B, 'FontWeight', 'Normal');
end

guidata(hObject, handles);
display_gonad(handles, hObject);
end

% --- Executes on button press in button_im3_R.
function button_im3_R_Callback(hObject, eventdata, handles)
% hObject    handle to button_im3_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if sum(handles.im3_index==1)
    handles.im3_index = setxor(handles.im3_index, 1);
    set(handles.button_im3_R, 'FontWeight', 'Normal');
    %     set(handles.button_im3_G, 'FontWeight', 'Normal');
    %     set(handles.button_im3_B, 'FontWeight', 'Normal');
else
    handles.im3_index = [handles.im3_index, 1];
    set(handles.button_im3_R, 'FontWeight', 'Bold');
    %     set(handles.button_im3_G, 'FontWeight', 'Normal');
    %     set(handles.button_im3_B, 'FontWeight', 'Normal');
end

guidata(hObject, handles);
display_gonad(handles, hObject);
end

% --- Executes on button press in button_im4_B.
function button_im4_B_Callback(hObject, eventdata, handles)
% hObject    handle to button_im4_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if sum(handles.im4_index==3)
    handles.im4_index = setxor(handles.im4_index, 3);
    %     set(handles.button_im4_R, 'FontWeight', 'Normal');
    %     set(handles.button_im4_G, 'FontWeight', 'Normal');
    set(handles.button_im4_B, 'FontWeight', 'Normal');
else
    handles.im4_index = [handles.im4_index, 3];
    %     set(handles.button_im4_R, 'FontWeight', 'Normal');
    %     set(handles.button_im4_G, 'FontWeight', 'Normal');
    set(handles.button_im4_B, 'FontWeight', 'Bold');
end

guidata(hObject, handles);
display_gonad(handles, hObject);

end

% --- Executes on button press in button_im4_G.
function button_im4_G_Callback(hObject, eventdata, handles)
% hObject    handle to button_im4_G (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if sum(handles.im4_index==2)
    handles.im4_index = setxor(handles.im4_index, 2);
    %     set(handles.button_im4_R, 'FontWeight', 'Normal');
    set(handles.button_im4_G, 'FontWeight', 'Normal');
    %     set(handles.button_im4_B, 'FontWeight', 'Normal');
else
    handles.im4_index = [handles.im4_index, 2];
    %     set(handles.button_im4_R, 'FontWeight', 'Normal');
    set(handles.button_im4_G, 'FontWeight', 'Bold');
    %     set(handles.button_im4_B, 'FontWeight', 'Normal');
end

guidata(hObject, handles);
display_gonad(handles, hObject);
end

% --- Executes on button press in button_im4_R.
function button_im4_R_Callback(hObject, eventdata, handles)
% hObject    handle to button_im4_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if sum(handles.im4_index==1)
    handles.im4_index = setxor(handles.im4_index, 1);
    set(handles.button_im4_R, 'FontWeight', 'Normal');
    %     set(handles.button_im4_G, 'FontWeight', 'Normal');
    %     set(handles.button_im4_B, 'FontWeight', 'Normal');
else
    handles.im4_index = [handles.im4_index, 1];
    set(handles.button_im4_R, 'FontWeight', 'Bold');
    %     set(handles.button_im4_G, 'FontWeight', 'Normal');
    %     set(handles.button_im4_B, 'FontWeight', 'Normal');
end

guidata(hObject, handles);
display_gonad(handles, hObject);

end




% ============================================================--
%
%  CHANNEL TEXTBOXES
%
% ============================================================--


function im1_scale_Callback(hObject, eventdata, handles)
% hObject    handle to im1_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of im1_scale as text
%        str2double(get(hObject,'String')) returns contents of im1_scale as a double

display_nucleus_vol_viz(handles);
end

% --- Executes during object creation, after setting all properties.
function im1_scale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to im1_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

function im1_off_Callback(hObject, eventdata, handles)
% hObject    handle to im1_off (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of im1_off as text
%        str2double(get(hObject,'String')) returns contents of im1_off as a double

display_nucleus_vol_viz(handles);
end

% --- Executes during object creation, after setting all properties.
function im1_off_CreateFcn(hObject, eventdata, handles)
% hObject    handle to im1_off (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function im1_gamma_Callback(hObject, eventdata, handles)
% hObject    handle to im1_gamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of im1_gamma as text
%        str2double(get(hObject,'String')) returns contents of im1_gamma as a double

display_nucleus_vol_viz(handles);
end

% --- Executes during object creation, after setting all properties.
function im1_gamma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to im1_gamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function im2_scale_Callback(hObject, eventdata, handles)
% hObject    handle to im2_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of im2_scale as text
%        str2double(get(hObject,'String')) returns contents of im2_scale as a double
display_nucleus_vol_viz(handles);
end

% --- Executes during object creation, after setting all properties.
function im2_scale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to im2_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function im2_off_Callback(hObject, eventdata, handles)
% hObject    handle to im2_off (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of im2_off as text
%        str2double(get(hObject,'String')) returns contents of im2_off as a double

display_nucleus_vol_viz(handles);
end

% --- Executes during object creation, after setting all properties.
function im2_off_CreateFcn(hObject, eventdata, handles)
% hObject    handle to im2_off (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function im2_gamma_Callback(hObject, eventdata, handles)
% hObject    handle to im2_gamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of im2_gamma as text
%        str2double(get(hObject,'String')) returns contents of im2_gamma as a double

display_nucleus_vol_viz(handles);
end

% --- Executes during object creation, after setting all properties.
function im2_gamma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to im2_gamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function im3_scale_Callback(hObject, eventdata, handles)
% hObject    handle to im3_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of im3_scale as text
%        str2double(get(hObject,'String')) returns contents of im3_scale as a double

display_nucleus_vol_viz(handles);
end

% --- Executes during object creation, after setting all properties.
function im3_scale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to im3_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function im3_off_Callback(hObject, eventdata, handles)
% hObject    handle to im3_off (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of im3_off as text
%        str2double(get(hObject,'String')) returns contents of im3_off as a double

display_nucleus_vol_viz(handles);
end

% --- Executes during object creation, after setting all properties.
function im3_off_CreateFcn(hObject, eventdata, handles)
% hObject    handle to im3_off (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function im3_gamma_Callback(hObject, eventdata, handles)
% hObject    handle to im3_gamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of im3_gamma as text
%        str2double(get(hObject,'String')) returns contents of im3_gamma as a double
display_nucleus_vol_viz(handles);

end

% --- Executes during object creation, after setting all properties.
function im3_gamma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to im3_gamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function im4_scale_Callback(hObject, eventdata, handles)
% hObject    handle to im4_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

display_nucleus_vol_viz(handles);

end
% --- Executes during object creation, after setting all properties.
function im4_scale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to im4_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

function im4_off_Callback(hObject, eventdata, handles)
% hObject    handle to im4_off (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
display_nucleus_vol_viz(handles);
end

% --- Executes during object creation, after setting all properties.
function im4_off_CreateFcn(hObject, eventdata, handles)
% hObject    handle to im4_off (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function im4_gamma_Callback(hObject, eventdata, handles)
% hObject    handle to im4_gamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
display_nucleus_vol_viz(handles);

end

% --- Executes during object creation, after setting all properties.
function im4_gamma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to im4_gamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% ============================================================--
%
%  Misc functions and object creation functions
%
% ============================================================--


function seg_thresh_Callback(hObject, eventdata, handles)
% hObject    handle to seg_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of seg_thresh as text
%        str2double(get(hObject,'String')) returns contents of seg_thresh as a double

end

% --- Executes during object creation, after setting all properties.
function seg_thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to seg_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function seg_smoothing_Callback(hObject, eventdata, handles)
% hObject    handle to seg_smoothing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of seg_smoothing as text
%        str2double(get(hObject,'String')) returns contents of seg_smoothing as a double

end

% --- Executes during object creation, after setting all properties.
function seg_smoothing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to seg_smoothing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in radio_intensity.
function radio_intensity_Callback(hObject, eventdata, handles)
% hObject    handle to radio_intensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_intensity
end

% --- Executes on button press in radio_orientation.
function radio_orientation_Callback(hObject, eventdata, handles)
% hObject    handle to radio_orientation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_orientation
end

% --- Executes on button press in radio_random.
function radio_random_Callback(hObject, eventdata, handles)
% hObject    handle to radio_random (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_random
end

function seg_channel_Callback(hObject, eventdata, handles)
% hObject    handle to seg_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of seg_channel as text
%        str2double(get(hObject,'String')) returns contents of seg_channel as a double
end

% --- Executes during object creation, after setting all properties.
function seg_channel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to seg_smoothing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes during object creation, after setting all properties.
function text55_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function initialize_channel_buttons(handles)

end


function listbox_surfs_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

% --- Executes during object creation, after setting all properties.
function alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


% --- Executes on button press in checkbox_use_smoothing.
function checkbox_use_smoothing_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_use_smoothing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_use_smoothing
end


function seg_gamma_Callback(hObject, eventdata, handles)
% hObject    handle to seg_gamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of seg_gamma as text
%        str2double(get(hObject,'String')) returns contents of seg_gamma as a double

end

% --- Executes during object creation, after setting all properties.
function seg_gamma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to seg_gamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in checkbox_displaySurfs.
function checkbox_displaySurfs_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_displaySurfs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_displaySurfs
end


% --- Executes on button press in checkbox_useMask.
function checkbox_useMask_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_useMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_useMask
end



% --- Executes during object creation, after setting all properties.
function textBox_dir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textBox_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


function ftrace_thresh_Callback(hObject, eventdata, handles)
% hObject    handle to ftrace_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ftrace_thresh as text
%        str2double(get(hObject,'String')) returns contents of ftrace_thresh as a double
end


% --- Executes during object creation, after setting all properties.
function ftrace_thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ftrace_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function ftrace_smoothing_Callback(hObject, eventdata, handles)
% hObject    handle to ftrace_smoothing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ftrace_smoothing as text
%        str2double(get(hObject,'String')) returns contents of ftrace_smoothing as a double
end


% --- Executes during object creation, after setting all properties.
function ftrace_smoothing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ftrace_smoothing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function ftrace_axis_channel_Callback(hObject, eventdata, handles)
% hObject    handle to ftrace_axis_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ftrace_axis_channel as text
%        str2double(get(hObject,'String')) returns contents of ftrace_axis_channel as a double
end


% --- Executes during object creation, after setting all properties.
function ftrace_axis_channel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ftrace_axis_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function ftrace_maxAxisLength_Callback(hObject, eventdata, handles)
% hObject    handle to ftrace_maxAxisLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ftrace_maxAxisLength as text
%        str2double(get(hObject,'String')) returns contents of ftrace_maxAxisLength as a double
end


% --- Executes during object creation, after setting all properties.
function ftrace_maxAxisLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ftrace_maxAxisLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function focus_thresh_Callback(hObject, eventdata, handles)
% hObject    handle to focus_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of focus_thresh as text
%        str2double(get(hObject,'String')) returns contents of focus_thresh as a double
end


% --- Executes during object creation, after setting all properties.
function focus_thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to focus_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function focus_smoothing_1_Callback(hObject, eventdata, handles)
% hObject    handle to focus_smoothing_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of focus_smoothing_1 as text
%        str2double(get(hObject,'String')) returns contents of focus_smoothing_1 as a double
end


% --- Executes during object creation, after setting all properties.
function focus_smoothing_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to focus_smoothing_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function focus_channel_1_Callback(hObject, eventdata, handles)
% hObject    handle to focus_channel_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of focus_channel_1 as text
%        str2double(get(hObject,'String')) returns contents of focus_channel_1 as a double
end



% --- Executes during object creation, after setting all properties.
function focus_channel_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to focus_channel_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end



function fTrace_focusChannel_Callback(hObject, eventdata, handles)
% hObject    handle to fTrace_focusChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fTrace_focusChannel as text
%        str2double(get(hObject,'String')) returns contents of fTrace_focusChannel as a double

end

% --- Executes during object creation, after setting all properties.
function fTrace_focusChannel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fTrace_focusChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in checkbox_displayFinalSeg.
function checkbox_displayFinalSeg_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_displayFinalSeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_displayFinalSeg
end



function focus_thresh_1_Callback(hObject, eventdata, handles)
% hObject    handle to focus_thresh_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of focus_thresh_1 as text
%        str2double(get(hObject,'String')) returns contents of focus_thresh_1 as a double

end

% --- Executes during object creation, after setting all properties.
function focus_thresh_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to focus_thresh_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function focus_smoothing_2_Callback(hObject, eventdata, handles)
% hObject    handle to focus_smoothing_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of focus_smoothing_2 as text
%        str2double(get(hObject,'String')) returns contents of focus_smoothing_2 as a double
end


% --- Executes during object creation, after setting all properties.
function focus_smoothing_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to focus_smoothing_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function focus_channel_2_Callback(hObject, eventdata, handles)
% hObject    handle to focus_channel_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of focus_channel_2 as text
%        str2double(get(hObject,'String')) returns contents of focus_channel_2 as a double
end

% --- Executes during object creation, after setting all properties.
function focus_channel_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to focus_channel_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function focus_thresh_2_Callback(hObject, eventdata, handles)
% hObject    handle to focus_thresh_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of focus_thresh_2 as text
%        str2double(get(hObject,'String')) returns contents of focus_thresh_2 as a double
end


% --- Executes during object creation, after setting all properties.
function focus_thresh_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to focus_thresh_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in checkbox_useSegGamma.
function checkbox_useSegGamma_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_useSegGamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_useSegGamma
end


% --- Executes on button press in checkbox_useOtsuThresh.
function checkbox_useOtsuThresh_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_useOtsuThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_useOtsuThresh
end


