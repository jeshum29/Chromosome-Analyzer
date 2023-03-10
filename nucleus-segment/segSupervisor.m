function varargout = segSupervisor(varargin)
% SEGSUPERVISOR MATLAB code for segSupervisor.fig
%      SEGSUPERVISOR, by itself, creates a new SEGSUPERVISOR or raises the existing
%      singleton*.
%
%      H = SEGSUPERVISOR returns the handle to a new SEGSUPERVISOR or the handle to
%      the existing singleton*.
%
%      SEGSUPERVISOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEGSUPERVISOR.M with the given input arguments.
%
%      SEGSUPERVISOR('Property','Value',...) creates a new SEGSUPERVISOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before segSupervisor_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to segSupervisor_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help segSupervisor

% Last Modified by GUIDE v2.5 22-Jul-2013 18:06:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @segSupervisor_OpeningFcn, ...
                   'gui_OutputFcn',  @segSupervisor_OutputFcn, ...
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

% --- Executes just before segSupervisor is made visible.
function segSupervisor_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to segSupervisor (see VARARGIN)

% Choose default command line output for segSupervisor

handles.output = hObject;

if length(varargin)==1
    handles.dd = varargin{1};
else
    handles.dd = get(handles.textBox_dir,'String');
end

handles.gonad = load([handles.dd filesep 'gonad.mat']);
handles.gonad = handles.gonad.gonad;
handles.showLabels = 1;

if ~isfield(handles.gonad, 'deleted_regions')
   handles.gonad.deleted_regions = [];
end

try
    handles.bandpass = load3([handles.gonad.writeDir filesep 'mk_bandpass.tif']);
catch
    makeBackgroundMask(handles.gonad);
    handles.bandpass = load3([handles.gonad.writeDir filesep 'mk_bandpass.tif']);
end


try
    handles.mkF = load3([handles.gonad.writeDir filesep 'mkF3_reseg.tif']);
    handles.pp = regionprops(handles.mkF,'Area','BoundingBox','Centroid');
    handles.mkFC = makeColorMask3(double(handles.mkF));
catch

    try
        handles.mkF = load3([handles.gonad.writeDir filesep 'mkF3.tif']);
        handles.pp = regionprops(handles.mkF,'Area','BoundingBox','Centroid');
        handles.mkFC = makeColorMask3(double(handles.mkF));
    catch
        handles.pp = [];
    end
    
end

% handles.im = load3([handles.gonad.writeDir filesep 'im' num2str(handles.gonad.SEG_CHANNEL) '.tif']);

handles.mkF_reseg = [];
handles.pp_reseg = [];

if ~isempty(handles.pp)
for ii = 1:length(handles.pp)
    handles.nuke_inds{ii} = num2str(ii);
end

set(handles.listbox_nukes, 'String', handles.nuke_inds, 'Value', 1);
end

set(handles.text_dir_name, 'String', handles.gonad.writeDir);


set(handles.edit_lowBandpassThresh,  'String', num2str(handles.gonad.BANDPASS_HIGH));
set(handles.edit_highBandpassThresh, 'String', num2str(handles.gonad.BANDPASS_LOW));
set(handles.edit_channel,            'String', num2str(handles.gonad.SEG_CHANNEL));

set(handles.edit_lowThresh, 'String', num2str(handles.gonad.LOW_THRESH/1000));
handles.gonad.HIGH_THRESH = handles.gonad.LOW_THRESH;

set(handles.edit_highThresh,  'String', num2str(handles.gonad.HIGH_THRESH/1000));
set(handles.edit_minVol,      'String', num2str(handles.gonad.MIN_AXIS_VOLUME));
set(handles.edit_maxVol,      'String', num2str(handles.gonad.MAX_AXIS_VOLUME));
set(handles.edit_maxDia,      'String', num2str(handles.gonad.MAX_NUKE_DIA));
set(handles.edit_aspectRatio, 'String', num2str(1));
 
set(handles.edit_reseg_lowThresh,   'String', num2str(handles.gonad.LOW_THRESH/1000));
set(handles.edit_reseg_highThresh,  'String', num2str(handles.gonad.HIGH_THRESH/1000));
set(handles.edit_reseg_minVol,      'String', num2str(handles.gonad.MIN_AXIS_VOLUME));
set(handles.edit_reseg_maxVol,      'String', num2str(handles.gonad.MAX_AXIS_VOLUME));
set(handles.edit_reseg_maxDia,      'String', num2str(handles.gonad.MAX_NUKE_DIA));
set(handles.edit_reseg_aspectRatio, 'String', num2str(1));

set(handles.edit_makeNukeData_minVol, 'String', num2str(handles.gonad.MIN_NUKE_VOLUME));
set(handles.edit_makeNukeData_maxVol, 'String', num2str(handles.gonad.MAX_NUKE_VOLUME));
set(handles.edit_makeNukeData_maxDia, 'String', num2str(handles.gonad.MAX_NUKE_DIA));

if ~isfield(handles.gonad, 'MAX_BORDER_AREA')
    handles.gonad.MAX_BORDER_AREA = 50;
end

set(handles.edit_makeNukeData_maxBorderArea, 'String', num2str(handles.gonad.MAX_BORDER_AREA));

if ~isfield(handles.gonad, 'NUCLEUS_CROP')
    handles.gonad.NUCLEUS_CROP = 5;
end

set(handles.edit_makeNukeData_crop, 'String', num2str(handles.gonad.NUCLEUS_CROP));

handles.index_selected = 1;

try
    display_image(handles,hObject);
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes segSupervisor wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = segSupervisor_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

end


function display_image(handles, hObject)

ff = figure(1);clf;
set(ff, 'WindowButtonDownFcn', @clicker);
set(ff, 'WindowKeyPressFcn', @presser);

set(gca,'Position',[0,0,1,1]);

screen_sz = getScreenSize;

if size(handles.mkF,1) > 512
    if screen_sz(3)==2560
        set(gcf, 'position', [1400 80 1140 1140]);
    else
        set(gcf, 'position', [900 0 1000 1000]);
    end
end

% set(gcf, 'position', [1402         373        1108        1109]);

display_image_internal(handles);
guidata(hObject, handles);

function display_image_internal(handles)

index_selected = handles.index_selected;

imb = squeeze(max(handles.mkFC,[],3));
mk = max(handles.mkF==index_selected,[],3);

imR = imb(:,:,1);
imG = imb(:,:,2);
imB = imb(:,:,3);

imR(mk) = 255/2;
imG(mk) = 255/2;
imB(mk) = 255/2;

imshow(cat(3, imR, imG, imB),[]);
hold on

FONT_SZ = 24;

if handles.showLabels
    
    MIN_NUKE_VOLUME  = str2double(get(handles.edit_makeNukeData_minVol,'String'));
    MAX_NUKE_VOLUME  = str2double(get(handles.edit_makeNukeData_maxVol,'String'));
    MAX_NUKE_DIA     = str2double(get(handles.edit_makeNukeData_maxDia,'String'));
    MAX_BORDER_AREA  = str2double(get(handles.edit_makeNukeData_maxBorderArea, 'String'));
    
    
    for ii = 1:length(handles.pp)
        
        if isempty(handles.pp(ii).Centroid)
            continue;
        end
        
        [rr,cc,zz] = getBoundingBox3(handles.pp, ii, 0, size(handles.mkF));
        
        dia = max([rr(end) - rr(1), cc(end) - cc(1), zz(end) - zz(1)]);
        
        [rr,cc,zz] = getBoundingBox3(handles.pp, ii, 3, size(handles.mkF));
        
        mkk = handles.mkF(rr,cc,zz)==ii;
        
        include_flag = isNuke2(...
            handles.pp(ii).Area, dia, mkk,...
            MIN_NUKE_VOLUME,...
            MAX_NUKE_VOLUME,...
            MAX_NUKE_DIA,...
            MAX_BORDER_AREA);

        if include_flag
            color = [50,255,50]/255;
            color = [1,1,1];
        else
            color = [200,10,10]/255;
        end
            
        if ~isempty(intersect(handles.gonad.deleted_regions, ii))
            color = [1,1,1]/2;
        end
            
        if ii==index_selected

            text(...
                handles.pp(ii).Centroid(1),...
                handles.pp(ii).Centroid(2),...
                num2str(ii),...
                'fontsize', round(FONT_SZ*1.5),'fontweight', 'bold', 'color', color);
        
            text(...
                20, 20,  ['Dia: ' num2str(dia) '   Area: ' num2str(handles.pp(ii).Area)],...
                'fontsize', FONT_SZ*1.5, 'fontweight', 'bold', 'color', [1,1,1]);
        else
            
            text(...
                handles.pp(ii).Centroid(1),...
                handles.pp(ii).Centroid(2),...
                num2str(ii), 'fontsize', FONT_SZ,'fontweight', 'bold', 'color', color);
        
        end
    end
    
end

end

function clicker(src, event)
    
    click_point = get(gca,'CurrentPoint');
    click_point = click_point(1,1:2);
    
    for ii = 1:length(handles.pp)
        
        if isempty(handles.pp(ii).Centroid)
            dist(ii) = 1e9;
            continue;
        end
        
        dist(ii) = sum((handles.pp(ii).Centroid(1:2) - click_point).^2);
        
    end
    
    minn = min(dist);
    ind = find(dist==minn);
    
    if ind==handles.index_selected
        return;
    end
    
    set(handles.listbox_nukes, 'Value', ind);
    handles.index_selected = ind;
    
    guidata(hObject, handles);
    display_image_internal(handles);
    display_nuke(handles);
    figure(1);
    
end

   
    function presser(src, event)
        
            if char(event.Character) == 's'
                handles.do_split = 0;
                button_split_Callback(hObject, [], handles);
            end
            
            if char(event.Character) == 'a'
                handles.do_split = 1;
                button_split_Callback(hObject, [], handles);
            end
            
            if char(event.Character) == 'd'
                button_delete_Callback(hObject, [], handles)
            end
        
    end

end



function handles = do_seg(handles)

handles.mkF   = makeNukesFromMask4(handles.bandpass, handles.gonad);
handles.pp    = regionprops(handles.mkF,'Area','BoundingBox','Centroid');
handles.mkFC  = makeColorMask3(handles.mkF);

for ii = 1:length(handles.pp)
    handles.nuke_inds{ii} = num2str(ii);
end

set(handles.listbox_nukes, 'String', handles.nuke_inds, 'Value', 1);

end

function handles = do_reseg(handles)

[rr,cc,zz] = getBoundingBox3(handles.pp, handles.index_selected, 3, size(handles.mkF));

handles.mkF_reseg = makeNukesFromMask3(handles.bandpass(rr,cc,zz), handles.gonad.RESEG);
handles.mkF_reseg = handles.mkF_reseg.*(handles.mkF(rr,cc,zz)==handles.index_selected);

handles.pp_reseg = regionprops(handles.mkF_reseg, 'Area', 'BoundingBox', 'Centroid');

end

function display_nuke(handles)

index_selected = handles.index_selected;
[rr,cc,zz] = getBoundingBox3(handles.pp, index_selected, 3, size(handles.mkF));
mkk = handles.mkF(rr,cc,zz);

vizMask(mkk==index_selected, [1,1,1], 99,1);
mkk(mkk==index_selected) = 0;

vizMask(mkk, [1,0,0], 99, 0);
end

function display_reseg(handles)

inds = 1:length(handles.pp_reseg);
cc = colormap(lines(length(inds)));

for ii = 1:length(inds)
    
    vizMask(handles.mkF_reseg==inds(ii), cc(ii,:), 999, ii==1);
    
end
end

% --- Executes on button press in button_bandpass.
function button_bandpass_Callback(hObject, eventdata, handles)

    handles.gonad.BANDPASS_HIGH  = str2double(get(handles.edit_lowBandpassThresh, 'String'));
    handles.gonad.BANDPASS_LOW   = str2double(get(handles.edit_highBandpassThresh, 'String'));

    handles.gonad.SEG_CHANNEL    = str2double(get(handles.edit_channel, 'String'));
    
    tic;makeBackgroundMask(handles.gonad);toc;
    handles.bandpass = load3([handles.gonad.writeDir filesep 'mk_bandpass.tif']);
    guidata(hObject, handles);
    
end


% --- Executes on button press in button_showLabels.
function button_showLabels_Callback(hObject, eventdata, handles)
% hObject    handle to button_showLabels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.showLabels = ~handles.showLabels;

display_image(handles,hObject);
guidata(hObject, handles);

end


% --- Executes on button press in button_segment.
function button_segment_Callback(hObject, eventdata, handles)

handles.gonad.LOW_THRESH   = str2double(get(handles.edit_lowThresh, 'String'))*1000;
handles.gonad.HIGH_THRESH  = str2double(get(handles.edit_highThresh, 'String'))*1000;
handles.gonad.HIGH_THRESH  = handles.gonad.LOW_THRESH;

handles.gonad.MIN_AXIS_VOLUME  = str2double(get(handles.edit_minVol, 'String'));
handles.gonad.MAX_AXIS_VOLUME  = str2double(get(handles.edit_maxVol, 'String'));
handles.gonad.MAX_NUKE_DIA     = str2double(get(handles.edit_maxDia, 'String'));
handles.gonad.ASPECT_RATIO     = str2double(get(handles.edit_aspectRatio, 'String'));
 
set(handles.edit_reseg_lowThresh, 'String', num2str(handles.gonad.LOW_THRESH/1000));
set(handles.edit_reseg_highThresh, 'String', num2str(handles.gonad.HIGH_THRESH/1000));

handles = do_seg(handles);
handles.pp = regionprops(handles.mkF,'Area','BoundingBox','Centroid');
guidata(hObject, handles);
display_image(handles,hObject);

end

% --- Executes on button press in button_reseg.
function button_reseg_Callback(hObject, eventdata, handles)
% hObject    handle to button_reseg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.gonad.RESEG.LOW_THRESH   = str2double(get(handles.edit_reseg_lowThresh, 'String'))*1000;
handles.gonad.RESEG.HIGH_THRESH  = str2double(get(handles.edit_reseg_highThresh, 'String'))*1000;
handles.gonad.RESEG.HIGH_THRESH  = handles.gonad.RESEG.LOW_THRESH;

handles.gonad.RESEG.MIN_AXIS_VOLUME  = str2double(get(handles.edit_reseg_minVol, 'String'));
handles.gonad.RESEG.MAX_AXIS_VOLUME  = str2double(get(handles.edit_reseg_maxVol, 'String'));
handles.gonad.RESEG.MAX_NUKE_DIA     = str2double(get(handles.edit_reseg_maxDia, 'String'));
handles.gonad.RESEG.MIN_NUKE_DIA     = handles.gonad.MIN_NUKE_DIA;
handles.gonad.RESEG.ASPECT_RATIO     = str2double(get(handles.edit_reseg_aspectRatio, 'String'));

handles = do_reseg(handles);

display_reseg(handles);

guidata(hObject, handles);
end


% --- Executes on button press in button_applyReseg.
function button_applyReseg_Callback(hObject, eventdata, handles)

if isempty(handles.mkF_reseg)
    return;
end

mkF_reseg = handles.mkF_reseg;

[rr,cc,zz] = getBoundingBox3(handles.pp, handles.index_selected, 3, size(handles.mkF));
mkFF = handles.mkF(rr,cc,zz);

for ii = 1:length(handles.pp_reseg)
    handles.pp_reseg(ii).Centroid = handles.pp_reseg(ii).Centroid  + [cc(1), rr(1), zz(1)];
    handles.pp_reseg(ii).BoundingBox(1:3) = handles.pp_reseg(ii).BoundingBox(1:3) + [cc(1), rr(1), zz(1)];
end

mkF_reseg(handles.mkF_reseg==1) = handles.index_selected;
handles.pp(handles.index_selected) = handles.pp_reseg(1);

if length(handles.pp_reseg)==1

    overlap = nonzeros(unique(mkFF(~~mkF_reseg)));
    other_inds = setxor(overlap, handles.index_selected);
    
    for ii = 1:length(other_inds)
        handles.pp(other_inds(ii)).Centroid = [];
        handles.pp(other_inds(ii)).BoundingBox = [];
    end
    
else

    maxx = length(handles.pp);
    
    for ii = 2:length(handles.pp_reseg)
        mkF_reseg(handles.mkF_reseg==ii) = maxx + ii - 1;
        handles.pp(maxx+ii-1) = handles.pp_reseg(ii);
    end

end

mkFF(mkFF==handles.index_selected) = 0;
mkFF(~~mkF_reseg) = mkF_reseg(~~mkF_reseg);

handles.mkF(rr,cc,zz) = mkFF;

guidata(hObject, handles);

display_image(handles,hObject);
update_index_listbox(handles);

end

% --- Executes on button press in button_merge.
function button_merge_Callback(hObject, eventdata, handles)

index_selected = handles.index_selected;
index_to_merge = get(handles.listbox_nukes,'Value');

[rr,cc,zz] = getBoundingBox3(handles.pp, [index_selected,index_to_merge], 3, size(handles.mkF));

mkFF = handles.mkF(rr,cc,zz);

mkFF(mkFF==index_to_merge) = index_selected;

pp = regionprops(mkFF,'Centroid','BoundingBox');

handles.pp(index_selected).Centroid          = pp(index_selected).Centroid + [cc(1), rr(1), zz(1)];
handles.pp(index_selected).BoundingBox(1:3)  = pp(index_selected).BoundingBox(1:3) + [cc(1), rr(1), zz(1)];
handles.pp(index_selected).BoundingBox(4:6)  = pp(index_selected).BoundingBox(4:6);

handles.pp(index_to_merge).Centroid     = [];
handles.pp(index_to_merge).BoundingBox  = [];

handles.mkF(rr,cc,zz) = mkFF;

guidata(hObject, handles);

display_image(handles,hObject);
update_index_listbox(handles);


end


% --- Executes on button press in button_delete.
function button_delete_Callback(hObject, eventdata, handles)

index_selected = handles.index_selected;

if ~isempty(intersect(handles.gonad.deleted_regions, index_selected))
    handles.gonad.deleted_regions = [handles.gonad.deleted_regions, index_selected];
else
    handles.gonad.deleted_regions = setxor(handles.gonad.deleted_regions, index_selected);
end

guidata(hObject, handles);
display_image(handles,hObject);


end




% --- Executes on button press in button_split.
function button_split_Callback(hObject, eventdata, handles)


try 
    handles.do_split;
catch
    handles.do_split = 0;
end

index_selected = handles.index_selected;
[rr,cc,zz] = getBoundingBox3(handles.pp, [index_selected], 3, size(handles.mkF));
mkFF = handles.mkF(rr,cc,zz)==index_selected;
mkL = bwlabeln(mkFF, 6);

pp = regionprops(mkL, 'Centroid','BoundingBox');
colors = colormap(lines(length(pp)));
nuke_cen = handles.pp(index_selected).Centroid;
    
    
% plot visualization
if ~handles.do_split

    for ii = 1:length(pp)

        vizMask(mkL==ii, colors(ii,:), 123, ii==1);

        offset = 3;

        surf_cen = pp(ii).Centroid;

        unit_vec = surf_cen - nuke_cen;
        unit_vec = unit_vec ./ sqrt(sum(unit_vec.*unit_vec));

        surf_cen = surf_cen + offset * unit_vec;

        text(...
            surf_cen(1), ...
            surf_cen(2), ...
            surf_cen(3), ...
            [num2str(ii)], 'fontsize', 24, 'color', colors(ii,:)*.7);

        subnuke_inds{ii} = num2str(ii);

    end

    set(handles.listbox_nukes,'String', subnuke_inds,'Value',1);

    handles.do_split = 1;
    figure(1);
else
    
    ind_to_split = get(handles.listbox_nukes,'Value');
        
    for ii = 1:length(pp)
        
        if isnan(pp(ii).Centroid(1))
            continue;
        end
        
        pp(ii).Centroid = pp(ii).Centroid  + [cc(1), rr(1), zz(1)];
        pp(ii).BoundingBox(1:3) = pp(ii).BoundingBox(1:3) + [cc(1), rr(1), zz(1)];
        
        if ii~=ind_to_split
            mkL(mkL==ii) = 999;
        else
            mkL(mkL==ii) = 1001;
        end
    end
    
    mkL(mkL==999) = 1;
    mkL(mkL==1001) = 2;

    pp = regionprops(mkL, 'Area', 'Centroid','BoundingBox');
    
    if length(pp)<2
        handles.do_split = 0;
        guidata(hObject, handles);
        'resetting split functions: try again'
        return;
    end
    
    for ii = 1:length(pp)
        
        pp(ii).Centroid = pp(ii).Centroid  + [cc(1), rr(1), zz(1)];
        pp(ii).BoundingBox(1:3) = pp(ii).BoundingBox(1:3) + [cc(1), rr(1), zz(1)];
        
    end
    
    mkL(mkL==1) = 999;
    mkL(mkL==2) = 1001;
    
    mkL(mkL==999) = index_selected;
    handles.pp(index_selected) = pp(1);

    maxx = length(handles.pp);

    mkL(mkL==1001) = maxx+1;
    handles.pp(maxx+1) = pp(2);

    mkFF = handles.mkF(rr,cc,zz);
    mkFF(~~mkL) = mkL(~~mkL);
    handles.mkF(rr,cc,zz) = mkFF;

    handles.do_split = 0;
    
    guidata(hObject, handles);
    display_image(handles,hObject);
    update_index_listbox(handles);
    
    
end

guidata(hObject, handles);
    
end

function update_index_listbox(handles)

for ii = 1:length(handles.pp)
    handles.nuke_inds{ii} = num2str(ii);
end

set(handles.listbox_nukes, 'String', handles.nuke_inds, 'Value', 1);
end


% --- Executes on button press in button_saveReseg.
function button_saveReseg_Callback(hObject, eventdata, handles)

'Saving segmentation'

save3(uint16(handles.mkF), [handles.gonad.writeDir filesep], 'mkF3_reseg', 0);
mkFC = makeColorMask3(double(handles.mkF));

mkFC_ = squeeze(max(mkFC,[],3));
imwrite(mkFC_, [handles.gonad.writeDir filesep 'mkF3_reseg_color_proj.tif']);

save3c((mkFC), [handles.gonad.writeDir filesep], 'mkF3_reseg_color', 0);


handles.gonad.MIN_NUKE_VOLUME  = str2double(get(handles.edit_makeNukeData_minVol,'String'));
handles.gonad.MAX_NUKE_VOLUME  = str2double(get(handles.edit_makeNukeData_maxVol,'String'));
handles.gonad.MAX_NUKE_DIA     = str2double(get(handles.edit_makeNukeData_maxDia,'String'));
handles.gonad.MAX_BORDER_AREA  = str2double(get(handles.edit_makeNukeData_maxBorderArea,'String'));

crop = str2double(get(handles.edit_makeNukeData_crop,'String'));
handles.gonad.NUCLEUS_CROP = crop;


gonad = handles.gonad;
save([handles.gonad.writeDir filesep 'gonad.mat'], 'gonad');

'Finished saving'

guidata(hObject, handles);

end

% --- Executes on button press in button_reseg_showLabels.
function button_reseg_showLabels_Callback(hObject, eventdata, handles)

end

% --- Executes on button press in button_showMakeNukeData.
function button_showMakeNukeData_Callback(hObject, eventdata, handles)

MIN_NUKE_VOLUME  = str2double(get(handles.edit_makeNukeData_minVol,'String'));
MAX_NUKE_VOLUME  = str2double(get(handles.edit_makeNukeData_maxVol,'String'));
MAX_NUKE_DIA     = str2double(get(handles.edit_makeNukeData_maxDia,'String'));
MAX_BORDER_AREA  = str2double(get(handles.edit_makeNukeData_maxBorderArea,'String'));

imb = squeeze(max(handles.mkFC,[],3));

imR = imb(:,:,1);
imG = imb(:,:,2);
imB = imb(:,:,3);

figure(66);
set(gca,'Position',[0,0,1,1]);
imshow(cat(3, imR, imG, imB),[]);
hold on

FONT_SZ = round(48 * size(imR,1) / 1000);

sz = size(handles.mkF);

props = regionprops(handles.mkF, 'Area', 'BoundingBox','Centroid');

for ii = 1:length(props)
    
    [rr,cc,zz] = getBoundingBox3(props, ii, 0, sz);
    dia = max([rr(end) - rr(1), cc(end) - cc(1), zz(end) - zz(1)]);
    
    
    [rr,cc,zz] = getBoundingBox3(props, ii, 3, sz);
    mkk = handles.mkF(rr,cc,zz)==ii;
    
    include_flag = isNuke2(...
        props(ii).Area, dia, mkk,...
        MIN_NUKE_VOLUME,...
        MAX_NUKE_VOLUME,...
        MAX_NUKE_DIA,...
        MAX_BORDER_AREA);
    
    ss = [num2str(ii) ' : ' num2str(dia) ' : ' num2str(handles.pp(ii).Area)];
    ss = num2str(ii);

    if include_flag
        
        text(...
            props(ii).Centroid(1),...
            props(ii).Centroid(2),...
            ss, 'fontsize', FONT_SZ,'fontweight', 'bold', 'color', [0,.8,0]);
        
        else
            
        text(...
            props(ii).Centroid(1),...
            props(ii).Centroid(2),...
            ss, 'fontsize', FONT_SZ,'fontweight', 'bold', 'color', [0,0,1]);
        
        end
    end
    

end


% --- Executes on button press in button_makeNukeData.
function button_makeNukeData_Callback(hObject, eventdata, handles)

handles.gonad.MIN_NUKE_VOLUME  = str2double(get(handles.edit_makeNukeData_minVol,'String'));
handles.gonad.MAX_NUKE_VOLUME  = str2double(get(handles.edit_makeNukeData_maxVol,'String'));
handles.gonad.MAX_NUKE_DIA     = str2double(get(handles.edit_makeNukeData_maxDia,'String'));
handles.gonad.MAX_BORDER_AREA  = str2double(get(handles.edit_makeNukeData_maxBorderArea,'String'));


crop = str2double(get(handles.edit_makeNukeData_crop,'String'));
handles.gonad.NUCLEUS_CROP = crop;

makeNukeData(handles.gonad, crop);

gonad = handles.gonad;
save([handles.gonad.writeDir filesep 'gonad.mat'], 'gonad');

guidata(hObject, handles);

end


% --- Executes on button press in button_save_raw_seg.
function button_save_raw_seg_Callback(hObject, eventdata, handles)

save3(uint16(handles.mkF), [handles.gonad.writeDir filesep], 'mkF3', 0);
save3(handles.mkFC, [handles.gonad.writeDir filesep], 'mkF3_color',0);
gonad = handles.gonad;
save([handles.gonad.writeDir filesep 'gonad.mat'], 'gonad');
guidata(hObject, handles);

end


% --- Executes on button press in button_showThresh.
function button_showThresh_Callback(hObject, eventdata, handles)

figure(11);

handles.gonad.LOW_THRESH   = str2double(get(handles.edit_lowThresh, 'String'))*1000;
handles.gonad.HIGH_THRESH  = str2double(get(handles.edit_highThresh, 'String'))*1000;
handles.gonad.HIGH_THRESH  = handles.gonad.LOW_THRESH;

handles.gonad.MIN_AXIS_VOLUME  = str2double(get(handles.edit_minVol, 'String'));
handles.gonad.MAX_AXIS_VOLUME  = str2double(get(handles.edit_maxVol, 'String'));
handles.gonad.MAX_NUKE_DIA     = str2double(get(handles.edit_maxDia, 'String'));

clf;
imshow(autogain(cat(3,...
    max(handles.bandpass > handles.gonad.LOW_THRESH,[],3),...
    max(handles.bandpass > handles.gonad.HIGH_THRESH,[],3),...
    0*handles.bandpass(:,:,1))), []);


gonad = handles.gonad;
save([handles.gonad.writeDir filesep 'gonad.mat'], 'gonad');


guidata(hObject, handles);

end


% --- Executes on button press in button_optiThresh.
function button_optiThresh_Callback(hObject, eventdata, handles)

bp = handles.bandpass;

thresh = graythresh(bp);

thresh = thresh*max(bp(:));

handles.gonad.LOW_THRESH = thresh;
handles.gonad.HIGH_THRESH = handles.gonad.LOW_THRESH;

set(handles.edit_lowThresh, 'String', num2str(round(thresh/1000)));
set(handles.edit_highThresh, 'String', num2str(round(thresh/1000)));

figure(11);
clf;
imshow(autogain(cat(3,...
    max(handles.bandpass > handles.gonad.LOW_THRESH,[],3),...
    max(handles.bandpass > handles.gonad.HIGH_THRESH,[],3),...
    0*handles.bandpass(:,:,1))), []);

guidata(hObject, handles);

end



% --- Executes on button press in button_reseg_showThresh.
function button_reseg_showThresh_Callback(hObject, eventdata, handles)

[rr,cc,zz] = getBoundingBox3(handles.pp, handles.index_selected, 3, size(handles.mkF));

imshow3ck(...
    handles.bandpass(rr,cc,zz) > handles.gonad.RESEG.LOW_THRESH,...
    handles.bandpass(rr,cc,zz) > handles.gonad.RESEG.HIGH_THRESH,...
    0, [], 3, 111);

guidata(hObject, handles);

end


% --- Executes on selection change in listbox_nukes.
function listbox_nukes_Callback(hObject, eventdata, handles)

% handles.index_selected = get(handles.listbox_nukes,'Value');
% 
% display_image(handles,hObject);
% display_nuke(handles);

guidata(hObject, handles);
end




% --- Executes on button press in button_showImage.
function button_showImage_Callback(hObject, eventdata, handles)

index_selected = handles.index_selected;
[rr,cc,zz] = getBoundingBox3(handles.pp, index_selected, 6, size(handles.mkF));
imm = autogain(handles.im(rr,cc,zz));
imm = imm - mean(imm(:));

vizIm(imm,999);


end





function edit_makeNukeData_maxVol_Callback(hObject, eventdata, handles)
display_image(handles,hObject);

end


function edit_makeNukeData_minVol_Callback(hObject, eventdata, handles)
display_image(handles,hObject);

end


function edit_makeNukeData_maxDia_Callback(hObject, eventdata, handles)

display_image(handles,hObject);

end











% --- Executes during object creation, after setting all properties.
function listbox_nukes_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

function edit_lowThresh_Callback(hObject, eventdata, handles)

end

% --- Executes during object creation, after setting all properties.
function edit_lowThresh_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit_highThresh_Callback(hObject, eventdata, handles)

end

% --- Executes during object creation, after setting all properties.
function edit_highThresh_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

function edit_minVol_Callback(hObject, eventdata, handles)

end

% --- Executes during object creation, after setting all properties.
function edit_minVol_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

function edit_maxVol_Callback(hObject, eventdata, handles)

end

% --- Executes during object creation, after setting all properties.
function edit_maxVol_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

function edit_maxDia_Callback(hObject, eventdata, handles)

end

% --- Executes during object creation, after setting all properties.
function edit_maxDia_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

function edit_reseg_lowThresh_Callback(hObject, eventdata, handles)

end

% --- Executes during object creation, after setting all properties.
function edit_reseg_lowThresh_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit_reseg_highThresh_Callback(hObject, eventdata, handles)

end

% --- Executes during object creation, after setting all properties.
function edit_reseg_highThresh_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit_reseg_minVol_Callback(hObject, eventdata, handles)

end

% --- Executes during object creation, after setting all properties.
function edit_reseg_minVol_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit_reseg_maxVol_Callback(hObject, eventdata, handles)

end

% --- Executes during object creation, after setting all properties.
function edit_reseg_maxVol_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit_reseg_maxDia_Callback(hObject, eventdata, handles)

end

% --- Executes during object creation, after setting all properties.
function edit_reseg_maxDia_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function edit_lowBandpassThresh_Callback(hObject, eventdata, handles)

end

% --- Executes during object creation, after setting all properties.
function edit_lowBandpassThresh_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit_highBandpassThresh_Callback(hObject, eventdata, handles)

end

% --- Executes during object creation, after setting all properties.
function edit_highBandpassThresh_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function edit_channel_Callback(hObject, eventdata, handles)
% hObject    handle to edit_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_channel as text
%        str2double(get(hObject,'String')) returns contents of edit_channel as a double

end
% --- Executes during object creation, after setting all properties.
function edit_channel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

% --- Executes during object creation, after setting all properties.
function edit_makeNukeData_minVol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_makeNukeData_minVol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end




% --- Executes during object creation, after setting all properties.
function edit_makeNukeData_maxVol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_makeNukeData_maxVol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function edit_makeNukeData_crop_Callback(hObject, eventdata, handles)
% hObject    handle to edit_makeNukeData_crop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_makeNukeData_crop as text
%        str2double(get(hObject,'String')) returns contents of edit_makeNukeData_crop as a double
end

% --- Executes during object creation, after setting all properties.
function edit_makeNukeData_crop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_makeNukeData_crop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function edit_aspectRatio_Callback(hObject, eventdata, handles)
% hObject    handle to edit_aspectRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_aspectRatio as text
%        str2double(get(hObject,'String')) returns contents of edit_aspectRatio as a double

end

% --- Executes during object creation, after setting all properties.
function edit_aspectRatio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_aspectRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


function edit_reseg_aspectRatio_Callback(hObject, eventdata, handles)
% hObject    handle to edit_reseg_aspectRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_reseg_aspectRatio as text
%        str2double(get(hObject,'String')) returns contents of edit_reseg_aspectRatio as a double

end

% --- Executes during object creation, after setting all properties.
function edit_reseg_aspectRatio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_reseg_aspectRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes during object creation, after setting all properties.
function edit_makeNukeData_maxDia_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_makeNukeData_maxDia (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function edit_makeNukeData_maxBorderArea_Callback(hObject, eventdata, handles)
% hObject    handle to edit_makeNukeData_maxBorderArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_makeNukeData_maxBorderArea as text
%        str2double(get(hObject,'String')) returns contents of edit_makeNukeData_maxBorderArea as a double

end

% --- Executes during object creation, after setting all properties.
function edit_makeNukeData_maxBorderArea_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_makeNukeData_maxBorderArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
