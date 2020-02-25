function varargout = FeedbackGUI(varargin)
% FEEDBACKGUI MATLAB code for FeedbackGUI.fig
%      FEEDBACKGUI, by itself, creates a new FEEDBACKGUI or raises the existing
%      singleton*.
%
%      H = FEEDBACKGUI returns the handle to a new FEEDBACKGUI or the handle to
%      the existing singleton*.
%
%      FEEDBACKGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FEEDBACKGUI.M with the given input arguments.
%
%      FEEDBACKGUI('Property','Value',...) creates a new FEEDBACKGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FeedbackGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FeedbackGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FeedbackGUI

% Last Modified by GUIDE v2.5 14-Sep-2018 15:26:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FeedbackGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @FeedbackGUI_OutputFcn, ...
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


% --- Executes just before FeedbackGUI is made visible.
function FeedbackGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FeedbackGUI (see VARARGIN)

global stop_signal;
global feedback;
global drift;
global adaptive;

stop_signal = false;
feedback = false;
drift = false;
adaptive = false; 


signal_handle = handles.signal_figure;
energy_handle = handles.energy_figure;
phonons_handle = handles.phonons_figure;
phonons_real_handle = handles.phonons_real_figure;
voltage_handle = handles.voltage_figure;
parameters_handle = handles.parameters;

kalman_test(3, signal_handle, energy_handle, phonons_handle, phonons_real_handle, voltage_handle, parameters_handle);

% Choose default command line output for FeedbackGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FeedbackGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FeedbackGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in signal.
function signal_Callback(hObject, eventdata, handles)
% hObject    handle to signal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global stop_signal;

signal_handle = handles.signal_figure;
energy_handle = handles.energy_figure;
phonons_handle = handles.phonons_figure;
phonons_real_handle = handles.phonons_real_figure;
voltage_handle = handles.voltage_figure;
parameters_handle = handles.parameters;

button_state = get(hObject,'Value');
if button_state == get(hObject,'Max') % just pressed
    cla reset;
    kalman_test(1e5, signal_handle, energy_handle, phonons_handle, phonons_real_handle, voltage_handle, parameters_handle);
    
elseif button_state == get(hObject,'Min')
    stop_signal = true;
end
    

% --- Executes on button press in activate_feedback.
function activate_feedback_Callback(hObject, eventdata, handles)
% hObject    handle to activate_feedback (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global feedback;
button_state = get(hObject,'Value');
if button_state == get(hObject,'Max') % just pressed
    feedback = true;
elseif button_state == get(hObject,'Min')
    feedback = false;
end

% --- Executes on button press in activate_drift.
function activate_drift_Callback(hObject, eventdata, handles)
% hObject    handle to activate_drift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global drift;
button_state = get(hObject,'Value');
if button_state == get(hObject,'Max') % just pressed
    drift = true;
elseif button_state == get(hObject,'Min')
    drift = false;
end


% --- Executes on button press in adaptive.
function adaptive_Callback(hObject, eventdata, handles)
% hObject    handle to adaptive (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global adaptive;
button_state = get(hObject,'Value');
if button_state == get(hObject,'Max') % just pressed
    adaptive = true;
elseif button_state == get(hObject,'Min')
    adaptive = false;
end
