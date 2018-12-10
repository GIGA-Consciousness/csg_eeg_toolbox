function varargout = csg_main(varargin)
%__________________________________________________________________________
%     ________  _____________ ___ 
%    / __ _/  |/ ___/_/_  __/|__ \
%   / /_ / /| |\__ \ \ / /   __/ /
%  / __// ___ |__/ / // / _ / __/ 
% /_/  /_/  |_|___/_//_/ (_)____/ 
%
% fMRI Artefact removal and Sleep Scoring Toolbox, FASST.2
% http://www.montefiore.ulg.ac.be/~phillips/FASST.html
%__________________________________________________________________________
%
% CSG_MAIN M-file for csg_main.fig
%      CSG_MAIN, by itself, creates a new CSG_MAIN or raises the existing
%      singleton.
%
%      H = CSG_MAIN returns the handle to a new CSG_MAIN or the handle to
%      the existing singleton*.
%
%      CSG_MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CSG_MAIN.M with the given input arguments.
%
%      CSG_MAIN('Property','Value',...) creates a new CSG_MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before crc_main_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to csg_main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%__________________________________________________________________
% Copyright (C) 2009 Cyclotron Research Centre

% Written by Y. Leclercq & C. Phillips, 2008.
% Cyclotron Research Centre, University of Liege, Belgium
% $Id$

% Edit the above text to modify the response to help csg_main

% Last Modified by GUIDE v2.5 21-Oct-2016 13:04:24

% Display ASCII Welcome
disp('                                                                              ');
disp('     ________  _____________ ___                __________ _____              ');
disp('    / ____/  |/ ___/_/_  __/|__ \              / _____/__// ___/              ');
disp('   / /_ / /| |\__ \ \ / /   __/ /     __      / /    \_ \/ /____              ');
disp('  / __// ___ |__/ / // / _ / __/     /__/    / /______/ / /__/_/              ');
disp(' /_/  /_/  |_|___/_//_/ (_)____/            /______/___/_____/                ');
disp('                                                                              ');
disp(' fMRI Artefact removal and Sleep Scoring Toolbox, FASST.2 - CSG               ');
disp(' http://www.montefiore.ulg.ac.be/~phillips/FASST.html                         ');
disp(' An SPM12-compatible toolbox.                                                  ');
fprintf('\n');

% Check if SPM is available, and maybe more one day...
ok = check_installation;
if ~ok
    beep
    fprintf('INSTALLATION PROBLEM!');
    return
end

% Add the fieldtrip toolbox from SPM, if necessary
if ~exist('ft_defaults','file')
    addpath(fullfile(spm('Dir'),'external','fieldtrip'));
end
ft_defaults;

% Check for Signal Processing Toolbox
persistent flag_TBX
if isempty(flag_TBX)
    flag_TBX = license('checkout','signal_toolbox');
    if ~flag_TBX
        pth = fullfile(spm_str_manip(mfilename('fullpath'),'h'),'SPTfunctions');
        addpath(pth)
        disp(['warning: using freely distributed equivalent to filtering functions ', ...
          'as Signal Processing Toolbox is not available.']);
    end
end


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @csg_main_OpeningFcn, ...
    'gui_OutputFcn',  @csg_main_OutputFcn, ...
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

% --- Executes just before csg_main is made visible.
function csg_main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to csg_main (see VARARGIN)


[A] = imread('LOGO_Simple.png','BackgroundColor',0.94*[1 1 1]);
image(A)
axis off

% Choose default command line output for csg_main
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes csg_main wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = csg_main_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in 'Display one file'.
function push_dis_main_Callback(hObject, eventdata, handles)
% hObject    handle to push_dis_main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%delete(handles.figure1)
csg_dis_selchan;

% --- Executes on button press in 'Compare multiple files'.
function push_dis_cmp_Callback(hObject, eventdata, handles)
% hObject    handle to push_dis_cmp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% flags.multcomp=1;
% setappdata(hObject,'multcomp',1);
flags.multcomp=1;
dis_selchan(flags);

% --- Executes on button press in 'Concatenate two files'.
function push_concatenate_Callback(hObject, eventdata, handles)
% hObject    handle to push_concatenate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
crc_concatenate;

% --- Executes on button press in 'Compute spectral power'.
function push_disfrqcomp_Callback(hObject, eventdata, handles)
% hObject    handle to push_disfrqcomp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dis_frqcomp;

% --- Executes on button press in 'Display spectral power'.
function push_freqplot_Callback(hObject, eventdata, handles)
% hObject    handle to push_freqplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dis_frq;

% --- Executes on button press in push_freqplotstat.
function push_freqplotstat_Callback(hObject, eventdata, handles)
% hObject    handle to push_freqplotstat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dis_frq;

% --- Executes on button press in 'Chunking tool'.
function push_chunk_Callback(hObject, eventdata, handles)
% hObject    handle to push_chunk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
csg_chunks;

% % --- Executes on button press in 'Artefact'.
% function push_artf_Callback(hObject, eventdata, handles)
% % hObject    handle to push_chunk (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% CSG_artifact;

% --- Executes on button press in 'channels definition'.
function push_chandef_Callback(hObject, eventdata, handles)
% hObject    handle to push_chunk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
csg_chandef;

% --- Executes on button press in 'Process data'.
function push_process_Callback(hObject, eventdata, handles)
% hObject    handle to push_chunk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
csg_process;

%% SUBFUNCTION

function ok = check_installation
% function to check installation state of toolbox,
% particullarly the SPM path setup

ok = true;

% Check SPM installation
if exist('spm.m','file')
    [SPMver, SPMrel] = spm('Ver');
    if ~(strcmpi(SPMver,'spm8') && str2double(SPMrel)>8.5) && ...
            ~strcmpi(SPMver,'spm12')
        beep
        fprintf('\nERROR:\n')
        fprintf('\tThe *latest* version of SPM8 or SPM12b should be installed on your computer,\n')
        fprintf('\tand be available on MATLABPATH!\n\n')
        ok = false;
    end
else
    beep
    fprintf('\nERROR:\n')
    fprintf('\tThe *latest* version of SPM8 should be installed on your computer,\n')
    fprintf('\tand be available on MATLABPATH!\n\n')
    ok = false;
end

return





