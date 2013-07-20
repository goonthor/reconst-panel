function varargout = reconst_panel(varargin)
% RECONST_PANEL MATLAB code for reconst_panel.fig
%      RECONST_PANEL, by itself, creates a new RECONST_PANEL or raises the existing
%      singleton*.
%
%      H = RECONST_PANEL returns the handle to a new RECONST_PANEL or the handle to
%      the existing singleton*.
%
%      RECONST_PANEL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RECONST_PANEL.M with the given input arguments.
%
%      RECONST_PANEL('Property','Value',...) creates a new RECONST_PANEL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before reconst_panel_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to reconst_panel_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help reconst_panel

% Last Modified by GUIDE v2.5 24-May-2013 11:45:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @reconst_panel_OpeningFcn, ...
    'gui_OutputFcn',  @reconst_panel_OutputFcn, ...
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

% --- Executes just before reconst_panel is made visible.
function reconst_panel_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to reconst_panel (see VARARGIN)

% path(path, './toolbox_diffc/toolbox_diffc');
% path(path, './toolbox_diffc/toolbox_diffc/toolbox');
% path('./libsvm', path);
path('./libsvm-weights-3.17/matlab',path);

% Initialize variables
% Used for frame-by-frame to indicate which set of 4 frames we are on.
handles.currpanel = 1;
% Flag indicating analysis needs to be run again as a result of changes
% user has made in the GUI.
handles.rebuildmodel = 1;
% We don't need to rebuild the model, but the new points have changed so we
% should query the model with the new set of points.
handles.needtoquery = 1;
% Sample points
handles.pnts = -1;
% Sample values
handles.vals = -1;
% New points
handles.newpnts = -1;
% New values
handles.guesses = -1;
% Number of dimensions in points.
handles.pntdim = -1;
% Number of samples
handles.numpts = -1;
% Number of dimenions in values
handles.valdim = -1;
% Number of values (should equal numpts)
handles.numvals = -1;
% Dimensionality of new points (should equal pntdim)
handles.newpntsdim = -1;
% Is the data valid?
handles.valid = 0;
% Flag for mesh data
handles.ismesh = 0;
% Points file name
handles.filename = 'No file selected';
% Integer for display coordinate system
handles.dispID = 1;
% Integer for sample coordinate system
handles.sampID = 1;
% Color scale data
handles.colorbounds = [-1 1];
% Flag for keogram new points
handles.isKeo = 0;
% Flag for slice plane new points
handles.isSlice = 0;
% SVR model
handles.model = 0;
% Parameter accompanying SVR model
handles.modelfactor = 0;
% Vector of weights
handles.weights = 0;
% Flag indicating user specified weights
handles.hasweights = 0;

% Choose default command line output for reconst_panel
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes reconst_panel wait for user response (see UIRESUME)
% uiwait(handles.figure1);
setappdata(0, 'hReconstPanel', gcf);

% --- Outputs from this function are returned to the command line.
function varargout = reconst_panel_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.guesses;

function displaybutton_Callback(hObject, eventdata, handles)

% First a variety of checks are run on the input data to check
% if it is formatted correctly.
if handles.pntdim < 0 || handles.valdim < 0
    disp(handles.pntdim);
    disp(handles.valdim);
    disp('No data currently entered.');
    return
end

if handles.pntdim > 4
    disp('The dimension of the points must be less than or equal to 4.')
    handles.valid = 0;
    return
end

if handles.numpts ~= handles.numvals
    disp('The number of points is not equal to the number of values.');
    handles.valid = 0;
    return
end

if handles.valdim > handles.pntdim
    disp('The dimension of the values is greater than the points.')
    handles.valid = 0;
    return
end

if handles.pntdim == 3 && handles.valdim ~= 1 && handles.valdim ~= 3
    disp('The dimension of the values must be either 1 or 3 for 3-D points.')
    handles.valid = 0;
    return
end

if handles.pntdim == 4 && handles.valdim ~= 1 && handles.valdim ~= 3
    disp('The dimension of the values must be either 1 or 3 for 4-D points.')
    handles.valid = 0;
    return
end

if handles.pntdim ~= handles.newpntsdim
    disp('Sample point dimension and new point dimension do not match.')
    disp(num2str(handles.pntdim))
    disp(num2str(handles.newpntsdim))
    handles.valid = 0;
    return
end

if handles.hasweights && numel(handles.weights) ~= numel(handles.vals(:,1))
    disp('The specified weights file has an incorrect number of entries.')
    handles.valid = 0;
    return
end

% If the input passes all tests it is probably valid.
handles.valid = 1;

% This flag is for SVR or TriScat
flag = get(handles.triradiobutton, 'Value')+1;

% This code runs the analysis. If only visual options were changed, it
% should not run again.
if handles.valid && (handles.rebuildmodel || handles.needstoquery)
    fsize = extractFeatureSize(handles);
    
    options.fsize = fsize;
    options.flag = flag;
    options.xval = get(handles.xvalcheckbox, 'value');
    options.isWin = get(handles.windowcheckbox, 'value');
    options.split = str2double(get(handles.splitedit,'String'));
    options.modelfactor = handles.modelfactor;
    options.rebuildmodel = handles.rebuildmodel;
    
    % Four options for newpoints format.
    % 1.
    if handles.isKeo
        handles.newpnts = generateKeogramPoints(handles);
        [handles.guesses fs handles.model handles.modelfactor] = ...
            performAnalysis(handles.pnts,handles.vals,handles.newpnts,...
            handles.weights, handles.model, options);
    % 2.
    elseif handles.isSlice
        [handles.newpnts handles.isSlice] = generateSlicePoints(handles);
        [handles.guesses fs handles.model handles.modelfactor] = ...
            performAnalysis(handles.pnts,handles.vals,handles.newpnts,...
            handles.weights, handles.model, options);
    % 3.
    elseif handles.ismesh
        % Generate a list of points from mesh form.
        pntmatrix = getMeshNewPoints(handles.pntdim, handles.newpnts);
        
        % Run analysis.
        [handles.guesses fs handles.model handles.modelfactor] = ...
            performAnalysis(handles.pnts,handles.vals,pntmatrix,...
            handles.weights, handles.model, options);
    % 4.
    else
        [handles.guesses fs handles.model handles.modelfactor] = ...
            performAnalysis(handles.pnts,handles.vals,handles.newpnts,...
            handles.weights, handles.model, options);
    end
    
    % Update feature size text.
    if numel(fs)>0
        set(handles.featuresizex1,'String', fs(1));
        if handles.pntdim >1
            set(handles.featuresizex2,'String', fs(2));
        end
        if handles.pntdim >2
            set(handles.featuresizex3,'String', fs(3));
        end
        if handles.pntdim >3
            set(handles.featuresizex4,'String', fs(4));
        end
    end
    
    handles.rebuildmodel = 0;
    handles.needstoquery = 0;
    guidata(hObject,handles);

else
    if ~handles.valid
        disp('Cannot update model: data is currently invalid');
    end
end

% If data is scalar, update color range.
if handles.valdim == 1
    colormin = get(handles.colorminedit, 'String');
    if numel(colormin) == 0
        handles.colorbounds(1) = min(handles.guesses);
        set(handles.colorminedit,'String', handles.colorbounds(1));
    else
        handles.colorbounds(1) = str2double(colormin);
    end
    colormax = get(handles.colormaxedit, 'String');
    if numel(colormax) == 0
        handles.colorbounds(2) = max(handles.guesses);
        set(handles.colormaxedit,'String', handles.colorbounds(2));
    else
        handles.colorbounds(2) = str2double(colormax);
    end
end

% Show plot.
if handles.valid
    visualize(handles.pnts, handles.vals, handles.newpnts, handles.guesses,...
        'showSamples', get(handles.showsampcheckbox,'Value'),...
        'showMag', get(handles.showmagcheckbox,'Value'),...
        'animate', get(handles.animateradiobutton, 'Value'),...
        'sbs', get(handles.sbsradiobutton, 'Value'),...
        'arrowsize', str2double(get(handles.arrowsizeedit,'String')),...
        'panel', handles.currpanel,...
        'fps', str2double(get(handles.fpsedit,'String')),...
        'dispID', handles.dispID,...
        'sampID', handles.sampID,...
        'title', handles.filename,...
        'keo', handles.isKeo,...
        'slice', handles.isSlice,...
        'cbounds', handles.colorbounds(1), handles.colorbounds(2));
else
    disp('Cannot visualize: data is currently invalid.');
end

guidata(hObject,handles);

function openptsbutton_Callback(hObject, eventdata, handles)
set(hObject, 'String', 'Opening...');
[filename,path,~] = uigetfile('*.mat');
set(hObject, 'String', filename);
disp(['Points file: ' filename])

pntsstruct = load([path filename]);
handles.filename = filename;
handles.pnts = pntsstruct.points;
handles.pntdim = size(handles.pnts,2);
handles.numpts = size(handles.pnts,1);
handles.rebuildmodel = 1;
handles.weights = 0;
handles.hasweights = 0;

disp(['The point dimension is ' num2str(handles.pntdim)]);
disp(['The number of points is ' num2str(handles.numpts)]);

if handles.pntdim >= 2
    set(handles.x2label,'enable','on');
    set(handles.featuresizex2,'enable','on');
    set(handles.polradiobutton,'enable','on');
    set(handles.polradiobuttonS,'enable','on');
    set(handles.animateradiobutton,'enable','on');
    set(handles.sbsradiobutton,'enable','on');
    set(handles.noanimradiobutton,'enable','on');
else
    set(handles.x2label,'enable','off');
    set(handles.featuresizex2,'enable','off');
    set(handles.featuresizex2,'String','');
    set(handles.polradiobutton,'enable','off');
    set(handles.polradiobuttonS,'enable','off');
    set(handles.animateradiobutton,'enable','off');
    set(handles.sbsradiobutton,'enable','off');
    set(handles.noanimradiobutton,'enable','off');
    set(handles.noanimradiobutton,'value',1);
    set(handles.fpsedit,'enable','off');
    set(handles.fpstext,'enable','off');
    set(handles.backbutton,'enable','off');
    set(handles.nextbutton,'enable','off');
end

if handles.pntdim >= 3
    set(handles.x3label,'enable','on');
    set(handles.featuresizex3,'enable','on');
    set(handles.sphradiobutton,'enable','on');
    set(handles.sphradiobuttonS,'enable','on');
    set(handles.keopushbutton, 'enable', 'on');
else
    set(handles.x3label,'enable','off');
    set(handles.featuresizex3,'enable','off');
    set(handles.featuresizex3,'String','');
    set(handles.sphradiobutton,'enable','off');
    set(handles.sphradiobuttonS,'enable','off');      
    set(handles.keopushbutton, 'enable', 'off');

end

if handles.pntdim == 4
    set(handles.x4label,'enable','on');
    set(handles.featuresizex4,'enable','on');
    set(handles.featuresizex2,'String','');
    set(handles.fpsedit,'enable','on');
    set(handles.fpstext,'enable','on');
    set(handles.animateradiobutton,'value',1);
    set(handles.noanimradiobutton,'enable','off');
    set(handles.triradiobutton,'enable','off');
    set(handles.svrradiobutton,'value',1);
    set(handles.slicepushbutton, 'enable', 'on');
else
    set(handles.x4label,'enable','off');
    set(handles.featuresizex4,'enable','off');
    set(handles.triradiobutton,'enable','on');
    set(handles.slicepushbutton, 'enable', 'off');
end

set(handles.cartradiobutton,'value',1);
set(handles.cartradiobuttonS,'value',1);

handles.dispID = 1;
handles.sampID = 1;
guidata(hObject,handles)

function openvalsbutton_Callback(hObject, eventdata, handles)
set(hObject, 'String', 'Opening...');
[filename,path,~] = uigetfile('*.mat');
set(hObject, 'String', filename);
disp(['Values file: ' filename])

valsstruct = load([path filename]);
handles.vals = valsstruct.values;
handles.valdim = size(handles.vals,2);
handles.numvals = size(handles.vals,1);
handles.rebuildmodel = 1;
handles.weights = 0;
handles.hasweights = 0;

disp(['The value dimension is ' num2str(handles.valdim)]);
disp(['The number of values is ' num2str(handles.numvals)]);

if handles.valdim == 2
    set(handles.irrotradiobutton,'enable','on');
    set(handles.incompradiobutton,'enable','on');
    set(handles.noneradiobutton,'enable','on');
    set(handles.showmagcheckbox, 'enable','on');
    set(handles.arrowtext, 'enable','on');
    set(handles.arrowsizeedit, 'enable','on');
else
    set(handles.irrotradiobutton,'enable','off');
    set(handles.incompradiobutton,'enable','off');
    set(handles.noneradiobutton,'enable','off');
    set(handles.noneradiobutton,'value',1);
    set(handles.showmagcheckbox, 'enable','off');
    set(handles.showmagcheckbox, 'value',0);
    set(handles.arrowtext, 'enable','off');
    set(handles.arrowsizeedit, 'enable','off');
end

guidata(hObject,handles)

function opennewbutton_Callback(hObject, eventdata, handles)
set(hObject, 'String', 'Opening...');
[filename,path,~] = uigetfile('*.mat');
set(hObject, 'String', filename);
disp(['New points file: ' filename])

pntsstruct = load([path filename]);
handles.needstoquery = 1;

if ismember('x1', fieldnames(pntsstruct))
    disp('New points in mesh form.');
    handles.newpnts = pntsstruct;
    handles.newpntsdim = numel(fieldnames(handles.newpnts));
    handles.ismesh = 1;
    switch handles.newpntsdim
        case 1
            handles.numpanels = ceil(numel(unique(handles.newpnts.x1))/4);
        case 2
            handles.numpanels = ceil(numel(unique(handles.newpnts.x2))/4);
        case 3
            handles.numpanels = ceil(numel(unique(handles.newpnts.x3))/4);
        case 4
            handles.numpanels = ceil(numel(unique(handles.newpnts.x4))/4);
    end
else
    disp('New points not in mesh form.')
    handles.newpnts = pntsstruct.newpoints;
    handles.newpntsdim = size(handles.newpnts,2);
    disp(['The number of new points is ' num2str(size(handles.newpnts,1))]);
    handles.ismesh = 0;
    switch handles.newpntsdim
        case 1
            handles.numpanels = ceil(numel(unique(handles.newpnts(:)))/4);
        case 2
            handles.numpanels = ceil(numel(unique(handles.newpnts(:,2)))/4);
        case 3
            handles.numpanels = ceil(numel(unique(handles.newpnts(:,3)))/4);
        case 4
            handles.numpanels = ceil(numel(unique(handles.newpnts(:,4)))/4);
    end 
end

handles.isKeo = 0;
handles.isSlice = 0;

set(handles.label1, 'Enable', 'off');
set(handles.label2, 'Enable', 'off');
set(handles.label1x1, 'Enable', 'off');
set(handles.editx1, 'Enable', 'off');
set(handles.label1x2, 'Enable', 'off');
set(handles.edity1, 'Enable', 'off');
set(handles.label2x1, 'Enable', 'off');
set(handles.editx2, 'Enable', 'off');
set(handles.label2x2, 'Enable', 'off');
set(handles.edity2, 'Enable', 'off');
set(handles.label1x3, 'Enable', 'off');
set(handles.editz1, 'Enable', 'off');
set(handles.label2x3, 'Enable', 'off');
set(handles.editz2, 'Enable', 'off');

disp(['The (new) point dimension is ' num2str(handles.newpntsdim)]);

guidata(hObject,handles)

function weightsbutton_Callback(hObject, eventdata, handles)
set(hObject, 'String', 'Opening...');
[filename,path,~] = uigetfile('*.mat');
set(hObject, 'String', filename);
disp(['Weights file: ' filename])

weightsstruct = load([path filename]);
handles.weights = weightsstruct.weights;
handles.hasweights = 1;
handles.rebuildmodel = 1;

guidata(hObject,handles)

function backbutton_Callback(hObject, eventdata, handles)
if handles.currpanel > 1
    handles.currpanel = handles.currpanel - 1;
end

if handles.valid
    visualize(handles.pnts, handles.vals, handles.newpnts, handles.guesses,...
        'showSamples', get(handles.showsampcheckbox,'Value'),...
        'showMag', get(handles.showmagcheckbox,'Value'),...
        'animate', get(handles.animateradiobutton, 'Value'),...
        'sbs', get(handles.sbsradiobutton, 'Value'),...
        'arrowsize', str2double(get(handles.arrowsizeedit,'String')),...
        'panel', handles.currpanel,...
        'fps', str2double(get(handles.fpsedit,'String')),...
        'dispID', handles.dispID,...
        'sampID', handles.sampID,...
        'title', handles.filename,...
        'cbounds', handles.colorbounds(1), handles.colorbounds(2));
else
    disp('Cannot visualize: data is currently invalid.');
end

guidata(hObject,handles);

function nextbutton_Callback(hObject, eventdata, handles)
if handles.currpanel < handles.numpanels
    handles.currpanel = handles.currpanel + 1;
end

if handles.valid
    visualize(handles.pnts, handles.vals, handles.newpnts, handles.guesses,...
        'showSamples', get(handles.showsampcheckbox,'Value'),...
        'showMag', get(handles.showmagcheckbox,'Value'),...
        'animate', get(handles.animateradiobutton, 'Value'),...
        'sbs', get(handles.sbsradiobutton, 'Value'),...
        'arrowsize', str2double(get(handles.arrowsizeedit,'String')),...
        'panel', handles.currpanel,...
        'fps', str2double(get(handles.fpsedit,'String')),...
        'dispID', handles.dispID,...
        'sampID', handles.sampID,...
        'title', handles.filename,...
        'cbounds', handles.colorbounds(1), handles.colorbounds(2));
else
    disp('Cannot visualize: data is currently invalid.');
end

guidata(hObject,handles);

function changealledit_Callback(hObject, eventdata, handles)
disp('Changing all the feature sizes.')
val = str2double(get(hObject,'String'));
if ~isnan(val)
    set(handles.featuresizex1,'String', val);
    if handles.pntdim >1
        set(handles.featuresizex2,'String', val);
    end
    if handles.pntdim >2
        set(handles.featuresizex3,'String', val);
    end
    if handles.pntdim >3
        set(handles.featuresizex4,'String', val);
    end
end
handles.rebuildmodel = 1;
guidata(hObject,handles);

function animatepanel_SelectionChangeFcn(hObject, eventdata, handles)
switch get(eventdata.NewValue,'Tag')
    case 'animateradiobutton'
        disp('Animate selected.')
        set(handles.nextbutton,'enable', 'off');
        set(handles.backbutton,'enable', 'off');
        set(handles.fpsedit,'enable', 'on');
        set(handles.fpstext,'enable', 'on');
    case 'sbsradiobutton'
        disp('Side-by-side selected.')
        set(handles.nextbutton,'enable', 'on');
        set(handles.backbutton,'enable', 'on');
        set(handles.fpsedit,'enable', 'off');
        set(handles.fpstext,'enable', 'off');
    case 'noanimradiobutton'
        disp('No animation.')
        set(handles.nextbutton,'enable', 'off');
        set(handles.backbutton,'enable', 'off');
        set(handles.fpsedit,'enable', 'off');
        set(handles.fpstext,'enable', 'off');
    otherwise
        disp('Something happened.')
end

function savebutton_Callback(hObject, eventdata, handles)
c = clock;
filename = ['results.' num2str(c(1)) '.' num2str(c(2)) '.' num2str(c(3))...
    '.' num2str(c(4)) '.' num2str(c(5)) '.mat'];
g = handles.guesses;
save(filename, 'g');

if handles.isSlice || handles.isKeo
    filename2 = ['newpoints.' num2str(c(1)) '.' num2str(c(2)) '.' num2str(c(3))...
    '.' num2str(c(4)) '.' num2str(c(5)) '.mat'];
    pnts = handles.newpnts;
    save(filename2, 'pnts');
end

if handles.ismesh
    filename2 = ['newpoints.' num2str(c(1)) '.' num2str(c(2)) '.' num2str(c(3))...
    '.' num2str(c(4)) '.' num2str(c(5)) '.mat'];
    pnts = getMeshNewPoints(handles.pntdim, handles.newpnts);
    save(filename2, 'pnts');
end

function resultsbutton_Callback(hObject, eventdata, handles)
[filename,path,~] = uigetfile('*.mat');
gstruct = load([path filename]);
handles.guesses = gstruct.g;
handles.rebuildmodel = 0;
handles.needstoquery = 0;
guidata(hObject,handles);

function coordpanel_SelectionChangeFcn(hObject, eventdata, handles)
switch get(eventdata.NewValue,'Tag')
    case 'cartradiobutton'
        handles.dispID = 1;
    case 'polradiobutton'
        handles.dispID = 2;
    case 'sphradiobutton'
        handles.dispID = 3;
    otherwise
        handles.dispID = 1;
end
guidata(hObject,handles);

function samplepanel_SelectionChangeFcn(hObject, eventdata, handles)
switch get(eventdata.NewValue,'Tag')
    case 'cartradiobuttonS'
        handles.sampID = 1;
    case 'polradiobuttonS'
        handles.sampID = 2;
    case 'sphradiobuttonS'
        handles.sampID = 3;
    otherwise
        handles.sampID = 1;
end
guidata(hObject,handles);

function keopushbutton_Callback(hObject, eventdata, handles)
if handles.valdim == 1
    handles.needstoquery = 1;
    handles.ismesh = 0;
    handles.newpntsdim = handles.pntdim;
    handles.isKeo = 1;
    handles.isSlice = 0;
    set(handles.label1, 'Enable', 'on');
    set(handles.label2, 'Enable', 'on');
    set(handles.label1x1, 'Enable', 'on');
    set(handles.editx1, 'Enable', 'on');
    set(handles.label1x2, 'Enable', 'on');
    set(handles.edity1, 'Enable', 'on');
    set(handles.label2x1, 'Enable', 'on');
    set(handles.editx2, 'Enable', 'on');
    set(handles.label2x2, 'Enable', 'on');
    set(handles.edity2, 'Enable', 'on');
    
    set(handles.keoslicepanel, 'Title', 'Keogram: specify end points of the segment');
    set(handles.label1, 'String', 'First Point');
    set(handles.label2, 'String', 'Second Point');
    
    if handles.pntdim == 4
        set(handles.label1x3, 'Enable', 'on');
        set(handles.editz1, 'Enable', 'on');
        set(handles.label2x3, 'Enable', 'on');
        set(handles.editz2, 'Enable', 'on');
    else
        set(handles.label1x3, 'Enable', 'off');
        set(handles.editz1, 'Enable', 'off');
        set(handles.label2x3, 'Enable', 'off');
        set(handles.editz2, 'Enable', 'off');
        
    end
else
    disp('Cannot show keogram with vector data.');
end
guidata(hObject,handles);

function slicepushbutton_Callback(hObject, eventdata, handles)
if handles.valdim == 1
    handles.needstoquery = 1;
    handles.ismesh = 0;
    handles.newpntsdim = handles.pntdim;
    handles.isKeo = 0;
    handles.isSlice = 1;
    set(handles.label1, 'Enable', 'on');
    set(handles.label2, 'Enable', 'on');
    set(handles.label1x1, 'Enable', 'on');
    set(handles.editx1, 'Enable', 'on');
    set(handles.label1x2, 'Enable', 'on');
    set(handles.edity1, 'Enable', 'on');
    set(handles.label2x1, 'Enable', 'on');
    set(handles.editx2, 'Enable', 'on');
    set(handles.label2x2, 'Enable', 'on');
    set(handles.edity2, 'Enable', 'on');
    set(handles.label1x3, 'Enable', 'on');
    set(handles.editz1, 'Enable', 'on');
    set(handles.label2x3, 'Enable', 'on');
    set(handles.editz2, 'Enable', 'on');
    
    set(handles.keoslicepanel, 'Title', 'Slice plane: specify point and normal.');
    set(handles.label1, 'String', 'Point');
    set(handles.label2, 'String', 'Normal'); 
else
    disp('Cannot show slice with vector data.');
end
guidata(hObject,handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following functions are for GUI elements which, if edited, 
% require the program to re-run the analysis.
function modelpanel_SelectionChangeFcn(hObject, eventdata, handles)
handles.rebuildmodel = 1;
guidata(hObject,handles);
function featuresizex1_Callback(hObject, eventdata, handles)
handles.rebuildmodel = 1;
guidata(hObject,handles);
function featuresizex2_Callback(hObject, eventdata, handles)
handles.rebuildmodel = 1;
guidata(hObject,handles);
function featuresizex3_Callback(hObject, eventdata, handles)
handles.rebuildmodel = 1;
guidata(hObject,handles);
function featuresizex4_Callback(hObject, eventdata, handles)
handles.rebuildmodel = 1;
guidata(hObject,handles);
function splitedit_Callback(hObject, eventdata, handles)
handles.rebuildmodel = 1;
guidata(hObject,handles);
function xvalcheckbox_Callback(hObject, eventdata, handles)
handles.rebuildmodel = 1;
guidata(hObject,handles);
function windowcheckbox_Callback(hObject, eventdata, handles)
handles.rebuildmodel = 1;
guidata(hObject,handles);
function vectorpanel_SelectionChangeFcn(hObject, eventdata, handles)
handles.rebuildmodel = 1;
guidata(hObject,handles);
function editx1_Callback(hObject, eventdata, handles)
handles.needstoquery = 1;
guidata(hObject,handles);
function edity1_Callback(hObject, eventdata, handles)
handles.needstoquery = 1;
guidata(hObject,handles);
function editz1_Callback(hObject, eventdata, handles)
handles.needstoquery = 1;
guidata(hObject,handles);
function editx2_Callback(hObject, eventdata, handles)
handles.needstoquery = 1;
guidata(hObject,handles);
function edity2_Callback(hObject, eventdata, handles)
handles.needstoquery = 1;
guidata(hObject,handles);
function editz2_Callback(hObject, eventdata, handles)
handles.needstoquery = 1;
guidata(hObject,handles);
%%%%%%%%%%%%%%%%%%%%%% END SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% UNUSED FUNCTIONS - Necessary for GUI. %%%%
%                  SKIP                       %
function showsampcheckbox_Callback(hObject, eventdata, handles)
function showmagcheckbox_Callback(hObject, eventdata, handles)
function featuresizex1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function featuresizex2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    disp('white')
    set(hObject,'BackgroundColor','white');
end
function featuresizex3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function featuresizex4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function changealledit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function arrowsizeedit_Callback(hObject, eventdata, handles)
function arrowsizeedit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function fpsedit_Callback(hObject, eventdata, handles)
function fpsedit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function splitedit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function colorminedit_Callback(hObject, eventdata, handles)
function colorminedit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function colormaxedit_Callback(hObject, eventdata, handles)
function colormaxedit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function editx1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edity1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function editz1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function editx2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edity2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function editz2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%%%%%%%%%%%%%%%%%%% END UNUSED FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%


