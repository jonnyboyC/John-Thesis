function varargout = mask_generation(varargin)
% MASK_GENERATION MATLAB code for mask_generation.fig
%      MASK_GENERATION, by itself, creates a new MASK_GENERATION or raises the existing
%      singleton*.
%
%      H = MASK_GENERATION returns the handle to a new MASK_GENERATION or the handle to
%      the existing singleton*.
%
%      MASK_GENERATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MASK_GENERATION.M with the given input arguments.
%
%      MASK_GENERATION('Property','Value',...) creates a new MASK_GENERATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mask_generation_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mask_generation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mask_generation

% Last Modified by GUIDE v2.5 13-Feb-2015 13:33:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mask_generation_OpeningFcn, ...
                   'gui_OutputFcn',  @mask_generation_OutputFcn, ...
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


% --- Executes just before mask_generation is made visible.
function mask_generation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mask_generation (see VARARGIN)

% Choose default command line output for mask_generation
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mask_generation wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mask_generation_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
