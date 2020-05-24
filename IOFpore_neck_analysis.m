function varargout = IOFpore_neck_analysis(varargin)
% IOFPORE_NECK_ANALYSIS MATLAB code for IOFpore_neck_analysis.fig
%      IOFPORE_NECK_ANALYSIS, by itself, creates a new IOFPORE_NECK_ANALYSIS or raises the existing
%      singleton*.
%
%      H = IOFPORE_NECK_ANALYSIS returns the handle to a new IOFPORE_NECK_ANALYSIS or the handle to
%      the existing singleton*.
%
%      IOFPORE_NECK_ANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IOFPORE_NECK_ANALYSIS.M with the given input arguments.
%
%      IOFPORE_NECK_ANALYSIS('Property','Value',...) creates a new IOFPORE_NECK_ANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before IOFpore_neck_analysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to IOFpore_neck_analysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help IOFpore_neck_analysis

% Last Modified by GUIDE v2.5 07-Mar-2018 21:06:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @IOFpore_neck_analysis_OpeningFcn, ...
                   'gui_OutputFcn',  @IOFpore_neck_analysis_OutputFcn, ...
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


% --- Executes just before IOFpore_neck_analysis is made visible.
function IOFpore_neck_analysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to IOFpore_neck_analysis (see VARARGIN)

% Choose default command line output for IOFpore_neck_analysis
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

selectFile(handles)

% UIWAIT makes IOFpore_neck_analysis wait for user response (see UIRESUME)
% uiwait(handles.figure1);


function selectFile(handles)
% Picks a new image file from a folder and displays the new image and
% Contrast/Brightness
[title, path] = uigetfile({strjoin({'*.tif','*.tiff','*.png','*.jpg'},';'),'All Image Files'});

global init
try init.Im = imread([path, title]);
    %init.Im = init.Im(:,:,:,1:3);
    set(handles.imageTitle,'String',[path,title])
catch
    init.Im = getimage(handles.axes);
    init.min = get(handles.CBmin,'Value');
    init.max = get(handles.CBmax,'Value');
    init.contrast = get(handles.CBcon,'Value');
    init.brightness = get(handles.CBbri,'Value');
    [path,title,ext] = fileparts(get(handles.imageTitle,'String'));
    path = [path,'/'];    title = [title,ext];
    set(handles.imageTitle,'String',[path,title])
    return
end
init.min = min(init.Im(:));
init.max = max(init.Im(:));
init.contrast = 0; %init.max - init.min;
init.brightness = mean(init.Im(:));

clear global results; hideResults(handles)
imshow(init.Im, 'Parent', handles.axes);
setValues(handles)


function setValues(handles)
global init
set(handles.CBmin,'Value',init.min);
set(handles.CBmax,'Value',init.max);
set(handles.CBcon,'Value',init.contrast);
set(handles.CBbri,'Value',init.brightness/2.55);


function updateImage(handles)
global init
% adjust lower threshold by minimum value
minI = get(handles.CBmin,'Value');
adjIm = init.Im; adjIm(adjIm<minI)=0;

% adjust upper threshold by maximum value 
maxI = get(handles.CBmax,'Value');
adjIm(adjIm>maxI)=255;%adjIm = imadjust(adjIm, [0 1], [0 maxI/255]);

% adjust contrast from -1 to 1 (starting at 0) and adjust in/out
% values by multiplying
contrast = get(handles.CBcon,'Value');
if contrast > 0
	adjIm = imadjust(adjIm, [0 1].*10^-contrast, [0 1]);
else
	adjIm = imadjust(adjIm, [0 1], [0 1].*10^contrast);
end

% adjust birghtness from 0 to 255 (starting at mean(Im(:)) ) and
% adjust in/out values to shift the new mean value by changing
% delta low-high
relBrightness = get(handles.CBbri,'Value')*2.55/mean(init.Im(:));
if relBrightness > 1
	adjIm = imadjust(adjIm, [0 1/relBrightness], [0 1]);
else
	adjIm = imadjust(adjIm, [1-relBrightness 1], [0 1]);
end

imshow(adjIm, 'Parent', handles.axes);
hideResults(handles);


% --- Executes on slider movement.
function CBmin_Callback(hObject, eventdata, handles)
% hObject    handle to CBmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(hObject,'Value');
set(hObject,'Value',round(value)); % round value to nearest integer
set(handles.CBminText,'String',num2str(round(value))); % update label
if get(handles.CBmax,'Value') <= value
    set(handles.CBmax,'Value',round(value+1)); % set CBmax to be larger than CBmin
end
updateImage(handles)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes on slider movement.
function CBmax_Callback(hObject, eventdata, handles)
% hObject    handle to CBmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(hObject,'Value');
set(hObject,'Value',round(value)); % round value to nearest integer
set(handles.CBmaxText,'String',num2str(round(value))); % update label
if get(handles.CBmin,'Value') >= value
    set(handles.CBmin,'Value',round(value-1)); % set CBmin to be smaller than CBmax
end
updateImage(handles)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes on slider movement.
function CBcon_Callback(hObject, eventdata, handles)
% hObject    handle to CBcon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.CBconText,'String',num2str(get(hObject,'Value'))); % update label
updateImage(handles)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes on slider movement.
function CBbri_Callback(hObject, eventdata, handles)
% hObject    handle to CBbri (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(hObject,'Value');
set(hObject,'Value',round(value)); % round value to nearest integer
set(handles.CBbriText,'String',num2str(round(value))); % update label
updateImage(handles)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider



% --- Executes on button press in checkPore.
function checkPore_Callback(hObject, eventdata, handles)
% hObject    handle to checkPore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.checkNeck,'Value',~get(hObject,'Value'))
set(handles.checkOpening,'Value',0)
updateImage(handles)
% Hint: get(hObject,'Value') returns toggle state of checkPore


% --- Executes on button press in checkOpening.
function checkOpening_Callback(hObject, eventdata, handles)
% hObject    handle to checkOpening (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.checkPore,'Value',0)
set(handles.checkNeck,'Value',~get(hObject,'Value'))
updateImage(handles)
% Hint: get(hObject,'Value') returns toggle state of checkOpening


% --- Executes on button press in checkNeck.
function checkNeck_Callback(hObject, eventdata, handles)
% hObject    handle to checkNeck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.checkPore,'Value',~get(hObject,'Value'))
set(handles.checkOpening,'Value',0)
updateImage(handles)
% Hint: get(hObject,'Value') returns toggle state of checkNeck


% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% Calculates the pore and neck size of close-packed inverse opal structure
% using Hough transform (for pores) and thresholding (for necks).
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global results
set(handles.runningText,'Visible','on'); drawnow;
if get(handles.checkPore,'Value')
Im = getimage(handles.axes);
    [centers,radii,metric] = imfindcircles(im2bw(Im),...
                [str2double(get(handles.poreMin,'String'))/2; str2double(get(handles.poreMax,'String'))/2], ...
                'Sensitivity', 0.95, ...
                'EdgeThreshold', 0.05, ...
                'Method', 'PhaseCode', ...
                'ObjectPolarity', 'dark');
    
	[X,Y] = meshgrid(1:numel(radii));
    overlap = any(reshape(sqrt((centers(X,1) - centers(Y,1)).^2 ...
                             + (centers(X,2) - centers(Y,2)).^2),size(X)) < radii(Y) ...
                  & metric(Y) > metric(X));
    centers(overlap,:) = []; radii(overlap) = []; metric(overlap) = [];
    
    results.poreCenters = centers(metric>0.2,:);
    
    if ~isempty(radii(metric>0.2))
        [IDX, D] = rangesearch(results.poreCenters,results.poreCenters, mean(radii(metric>0.2))*3);
        conn = arrayfun(@(n)arrayfun(@(x){[centers(n,1) centers(x,1)]; [centers(n,2) centers(x,2)]},IDX{n},'uni',0),1:length(IDX),'uni',0); conn = [conn{:}]; conn = [conn{:}];

        results.poreRadii = cellfun(@(d)mean(d(d>0)),D)/2;

        viscircles(handles.axes,results.poreCenters,results.poreRadii);
        hold(gca,'on'); plot(conn{:},'Color','r','LineWidth',2); hold(gca,'off')
    %     [~,D] = knnsearch(centers,centers,'k',7);
    %     dist = D(:); dist(dist == 0 | dist>1.5*median(dist)) = [];
        results.poreSize = nanmean(results.poreRadii)*2;
    end
elseif get(handles.checkOpening,'Value')
    Im = getimage(handles.axes);
    [centers,radii,metric] = imfindcircles(im2bw(Im),...
                [str2double(get(handles.poreMin,'String'))/2; str2double(get(handles.poreMax,'String'))/2], ...
                'Sensitivity', 0.95, ...
                'EdgeThreshold', 0.05, ...
                'Method', 'PhaseCode', ...
                'ObjectPolarity', 'dark');
    
	[X,Y] = meshgrid(1:numel(radii));
    overlap = any(reshape(sqrt((centers(X,1) - centers(Y,1)).^2 ...
                             + (centers(X,2) - centers(Y,2)).^2),size(X)) < radii(Y) ...
                  & metric(Y) > metric(X));
    centers(overlap,:) = []; radii(overlap) = []; metric(overlap) = [];
    
    results.openingCenters = centers(metric>0.2,:);
    results.openingRadii = radii(metric>0.2);
    
    viscircles(handles.axes,results.openingCenters,results.openingRadii);
    results.openingSize = mean(results.openingRadii);
elseif get(handles.checkNeck,'Value')
    Im = getimage(handles.axes);
    necks = regionprops(imcomplement(imbinarize(Im,0.1)), 'Orientation', 'MajorAxisLength', ...
    'MinorAxisLength', 'Eccentricity', 'Centroid');
    necks([necks.MajorAxisLength] < 10) = []; %removing isolated pixels (15 before)
    necks([necks.MajorAxisLength] < 2/3*median([necks.MajorAxisLength])) = []; %15 before
    necks([necks.MajorAxisLength] > 2*median([necks.MajorAxisLength])) = [];
    results.necks = necks;
    plotNecks(handles,necks)
    results.neckSize = mean([necks.MajorAxisLength]);
end
set(handles.runningText,'Visible','off');

% --- Executes on button press in show.
function show_Callback(hObject, eventdata, handles)
% Plots the pores or necks found as an overlay to the image.
% hObject    handle to show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global results

if ~isfield(results,'scale')
    set(handles.resultsText,'String','scale bar undefined');
    set(handles.measureScale,'BackgroundColor','r');
    return
end

set(handles.resultsText,'String','')

possibleResults = {'poreSize','neckSize','openingSize'};
found = ismember(possibleResults,fields(results));

if isempty(possibleResults(found))
    set(handles.resultsText,'String','no results yet')
else
    if get(handles.checkPore,'Value')
        viscircles(handles.axes,results.poreCenters,results.poreRadii);
    elseif get(handles.checkNeck,'Value')
        plotNecks(handles,results.necks)
    elseif get(handles.checkOpening,'Value')
        viscircles(handles.axes,results.openingCenters,results.openingRadii);
    end
    
    res = nan(size(possibleResults));
    for result = possibleResults(found)
        res(strcmp(possibleResults,result{:})) = round(results.(result{:})*results.scale);
    end
    set(handles.resultsText,'String', ...
        sprintf('Pore diameter: %3.0fnm\nNeck diameter: %3.1fnm\nOpening diameter: %3.0fnm',res));
end


% --- Executes on button press in export.
function export_Callback(hObject, eventdata, handles)
% exports the measured pore diameter and neck angle into a text file
% corresponding to the name of the image
% hObject    handle to export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global results

if ~isfield(results,'scale')
    set(handles.resultsText,'String','scale bar undefined');
    set(handles.measureScale,'BackgroundColor','r');
    return
end

possibleResults = {'poreSize','neckSize','openingSize'};
found = ismember(possibleResults,fields(results));

res = nan(size(possibleResults));
for result = possibleResults(found)
    res(strcmp(possibleResults,result{:})) = results.(result{:})*results.scale;
end

% switch isfield(results,'poreSize') + isfield(results,'neckSize')*2
%     case 1
%         res = [results.poreSize*results.scale];
%     case 2
%         res = [NaN results.neckSize*results.scale];
%     case 3
%         res = [results.poreSize*results.scale results.neckSize*results.scale];
%     otherwise
%         res = [];
% end

[path,title] = fileparts(get(handles.imageTitle,'String'));
filename = [path,'/',title,'.txt'];

fileID = fopen(filename,'w');
fprintf(fileID,'Pore diameter: %3.0fnm\nNeck diameter: %3.1fnm\nOpening diameter: %3.0fnm',res);
fclose(fileID);


% --- Executes on button press in adjustCB.
function adjustCB_Callback(hObject, eventdata, handles)
% hObject    handle to adjustCB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imcontrast()


% --- Executes on button press in changeImage.
function changeImage_Callback(hObject, eventdata, handles)
% hObject    handle to changeImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selectFile(handles)


% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global init
imshow(init.Im,'Parent',handles.axes)
setValues(handles)


% --- Executes on button press in measureScale.
function measureScale_Callback(hObject, eventdata, handles)
% hObject    handle to measureScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try [x,y,~] = improfile;
catch
    datacursormode toggle
    return
end
set(handles.scale,'String',num2str(max([abs(x(end)-x(1)), abs(y(end)-y(1))])));
set(handles.measureScale,'BackgroundColor',0.94*[1 1 1]);
set(handles.resultsText,'String','');

global results
results.scale = str2double(get(handles.scaleEdit,'String'))/str2double(get(handles.scale,'String'));


function scaleEdit_Callback(hObject, eventdata, handles)
% hObject    handle to scaleEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try str2double(get(hObject,'String'));
catch
    set(hObject,'String','200');
end

global results
results.scale = str2double(get(handles.scaleEdit,'String'))/str2double(get(handles.scale,'String'));
% Hints: get(hObject,'String') returns contents of scaleEdit as text
%        str2double(get(hObject,'String')) returns contents of scaleEdit as a double


% --- Executes on button press in measurePore.
function measurePore_Callback(hObject, eventdata, handles)
% hObject    handle to measurePore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try [x,y,~] = improfile;
catch
    datacursormode toggle
    return
end
radius = sqrt((x(end)-x(1))^2 + (y(end)-y(1))^2)/2;
min = round(radius-0.25*radius);
set(handles.poreMin,'String',num2str(min*2));
max = round(radius+0.25*radius);
set(handles.poreMax,'String',num2str(max*2));


function poreMin_Callback(hObject, eventdata, handles)
% hObject    handle to poreMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try str2double(get(hObject,'String'));
catch
    set(hObject,'String','60');
end
% Hints: get(hObject,'String') returns contents of poreMin as text
%        str2double(get(hObject,'String')) returns contents of poreMin as a double


function poreMax_Callback(hObject, eventdata, handles)
% hObject    handle to poreMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try str2double(get(hObject,'String'));
catch
    set(hObject,'String','100');
end
% Hints: get(hObject,'String') returns contents of poreMax as text
%        str2double(get(hObject,'String')) returns contents of poreMax as a double


%%% Unchanged functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Outputs from this function are returned to the command line.
function varargout = IOFpore_neck_analysis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function CBmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CBmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object creation, after setting all properties.
function CBmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CBmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object creation, after setting all properties.
function CBcon_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CBcon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object creation, after setting all properties.
function CBbri_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CBbri (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object creation, after setting all properties.
function scaleEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scaleEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function poreMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to poreMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function poreMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to poreMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%% Misc functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotNecks(handles,s)
% plotNecks(handles,s) displays the locations of necks over the image
% s    stats for the locations of neck obtained using regionprops 
%      (e.g. necks = regionprops(imcomplement(imbinarize(Im,0.1)), 'Orientation', 'MajorAxisLength', ...
%      'MinorAxisLength', 'Eccentricity', 'Centroid');
hold(handles.axes,'on')
phi = linspace(0,2*pi,50);
cosphi = cos(phi);
sinphi = sin(phi);


for k = 1:length(s)
    xbar = s(k).Centroid(1);
    ybar = s(k).Centroid(2);

    a = s(k).MajorAxisLength/2;
    b = s(k).MinorAxisLength/2;

    theta = pi*s(k).Orientation/180;
    R = [ cos(theta)   sin(theta)
         -sin(theta)   cos(theta)];

    xy = [a*cosphi; b*sinphi];
    xy = R*xy;

    x = xy(1,:) + xbar;
    y = xy(2,:) + ybar;

    plot(handles.axes,x,y,'r','LineWidth',2);
end
hold off
%set(handles.figure1,'Pointer','arrow')


function hideResults(handles)
% sets the results text empty.
set(handles.resultsText,'String','')
