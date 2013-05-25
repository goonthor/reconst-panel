function varargout = svrtester(varargin)
% SVRTESTER MATLAB code for svrtester.fig
%      SVRTESTER, by itself, creates a new SVRTESTER or raises the existing
%      singleton*.
%
%      H = SVRTESTER returns the handle to a new SVRTESTER or the handle to
%      the existing singleton*.
%
%      SVRTESTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SVRTESTER.M with the given input arguments.
%
%      SVRTESTER('Property','Value',...) creates a new SVRTESTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before svrtester_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to svrtester_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help svrtester

% Last Modified by GUIDE v2.5 21-Feb-2013 12:45:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @svrtester_OpeningFcn, ...
                   'gui_OutputFcn',  @svrtester_OutputFcn, ...
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


function svrtester_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
path('./libsvm', path);
guidata(hObject, handles);

function varargout = svrtester_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function ampedit_Callback(hObject, eventdata, handles)

function ampedit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function freqedit_Callback(hObject, eventdata, handles)

function freqedit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sampedit_Callback(hObject, eventdata, handles)

function sampedit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function gammaedit_Callback(hObject, eventdata, handles)

function gammaedit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function costedit_Callback(hObject, eventdata, handles)

function costedit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function nuedit_Callback(hObject, eventdata, handles)

function nuedit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function epedit_Callback(hObject, eventdata, handles)

function epedit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function guessedit_Callback(hObject, eventdata, handles)

function guessedit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% **** BEGIN RELEVANT CODE **** %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gobutton_Callback(hObject, eventdata, handles)
MIN = -1;
MAX = 1;
sampnum = str2double(get(handles.sampedit,'String'));
isPer = get(handles.perradiobutton, 'Value');

amp = str2double(get(handles.ampedit,'String'));
f = str2double(get(handles.freqedit,'String'));

newnum = str2double(get(handles.guessedit,'String'));

nu = str2double(get(handles.nuedit,'String'));
ep = str2double(get(handles.epedit,'String'));
g = str2double(get(handles.gammaedit,'String'));
c = str2double(get(handles.costedit,'String'));
isNu = get(handles.nuradiobutton, 'Value');
type = 3;
if isNu
    type = 4;
end

fnhandle = @(z) amp*sin(2*pi*f*z);

x = getSamples(sampnum, isPer, MIN, MAX);
y = fnhandle(x);
newpoints = linspace(MIN, MAX, newnum)';
vals = doSVR(x, y, newpoints, type, nu, ep, g, c);
plotorig(MIN, MAX, f, fnhandle);
graphresults(x, y, newpoints, vals);


function x = getSamples(sampnum, isPer, min, max)
if isPer
    x = linspace(min, max, sampnum)';
else
    x = (max-min)*rand(sampnum, 1) + min;
end
   
function plotorig(min, max, f, fnhandle)
x_orig = linspace(min, max, f*100);
figure(2)
clf(2)
hold on
title('Original Signal')
plot(x_orig, fnhandle(x_orig))
hold off

function vals = doSVR(x, y, newpoints, type, nu, ep, g, c)
options = ['-q -s ' num2str(type) ' -t 2 -g ' num2str(g) ' -c ' num2str(c) ' -n ' num2str(nu) ' -p ' num2str(ep)];
tic
model = svmtrain(y, x, options);
vals = svmpredict(rand(size(newpoints,1),1), newpoints, model);
t = toc;
disp(['Analysis completed in ' num2str(t) ' seconds.']);


function graphresults(x, y, newpoints, vals)
figure(1)
clf(1)
hold on
title('Reconstructed Signal (samples in red)')
plot(newpoints, vals)
plot(x, y, 'r*')
hold off


