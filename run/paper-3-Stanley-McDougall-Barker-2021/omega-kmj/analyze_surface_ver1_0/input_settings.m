function varargout = input_settings(varargin)

%           input_settings.m -- GUI for setting input settings.
%
% Usage:    analyse_surface.m
%
%           A graphical user interface (GUI) used to set input settings, 
%           such as Arakawa gird and wrap. Called from analyse_surface.m
%
%   _________________________________________________________________
%   This is part of the analyze_surface toolbox, (C) 2009 A. Klocker
%   type 'help analyze_surface' for more information 
%   type 'analyze_surface_license' for license details
%   type 'analyze_surface_version' for version details
%

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @input_settings_OpeningFcn, ...
                   'gui_OutputFcn',  @input_settings_OutputFcn, ...
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

function input_settings_OpeningFcn(hObject, eventdata, handles, varargin) %#ok

handles.output = hObject;
handles.ocean.grid = 'bgrid';
handles.ocean.wrap = 'none';

% Update handles structure

guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = input_settings_OutputFcn(hObject, eventdata, handles) %#ok

% Get default command line output from handles structure
uiwait

varargout{1} = handles.output;

%ocean_options = handles.ocean;

close(gcf)

%% Arakawa Grid

function bgrid_Callback(hObject, eventdata, handles) %#ok

handles.ocean.grid = 'bgrid';
guidata(hObject, handles);

%--------------------------------------------------------------------------

function cgrid_Callback(hObject, eventdata, handles) %#ok

handles.ocean.grid = 'cgrid';
guidata(hObject, handles);

%% Wrap

function wrap_none_Callback(hObject, eventdata, handles) %#ok

handles.ocean.wrap = 'none';
guidata(hObject, handles);

%--------------------------------------------------------------------------

function wrap_long_Callback(hObject, eventdata, handles) %#ok

handles.ocean.wrap = 'long';
guidata(hObject, handles);

%% OK

function go_ok_Callback(hObject, eventdata, handles) %#ok

handles.ocean

global ocean_options

ocean_options = handles.ocean;

uiresume(gcf)