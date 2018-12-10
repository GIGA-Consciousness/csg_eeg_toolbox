function varargout = csg_process(varargin)
% CSG_PROCESS MATLAB code for csg_process.fig
%      CSG_PROCESS, by itself, creates a new CSG_PROCESS or raises the existing
%      singleton*.
%
%      H = CSG_PROCESS returns the handle to a new CSG_PROCESS or the handle to
%      the existing singleton*.
%
%      CSG_PROCESS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CSG_PROCESS.M with the given input arguments.
%
%      CSG_PROCESS('Property','Value',...) creates a new CSG_PROCESS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before csg_process_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to csg_process_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help csg_process

% Last Modified by GUIDE v2.5 08-Nov-2016 11:08:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @csg_process_OpeningFcn, ...
                   'gui_OutputFcn',  @csg_process_OutputFcn, ...
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


% --- Executes just before csg_process is made visible.
function csg_process_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to csg_process (see VARARGIN)

% Choose default command line output for csg_process
handles.output = hObject;

prefile = spm_select(1, 'any', 'Select imported EEG file','' ...
    ,pwd,'\.[mMvVeErR][dDhHaA][fFDdTtwW]');
handles.multcomp=0;
D{1} = crc_eeg_load(deblank(prefile));
file = fullfile(D{1}.path,D{1}.fname);
handles.file = file;
handles.chan = upper(chanlabels(D{1}));
handles.Dmeg = D;

% chanselection for the power spectrum computation
%handles.chansel = {'All ''good'' EEG', 'Load Selection', 'Select channels'};
handles.epoch = 4;
handles.powchan = find(strcmp(upper(chantype(handles.Dmeg{1})),'EEG'));
% visibility 
set(handles.sel_badchan,'enable','on')
set(handles.sel_badepoch,'enable','on')
handles.flags_badchan = 0;
handles.flags_badepoch = 0;
handles.flags_interpol = 0;
handles.flags_reref = 0;
handles.flags_ps = 0;
handles.avgmas = 1;
% set(handles.list_selchan,'String',handles.chansel, 'Value', 1);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes csg_process wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = csg_process_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in sel_badchan.
function sel_badchan_Callback(hObject, eventdata, handles)
% hObject    handle to sel_badchan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sel_badchan
handles.flags_badchan = get(hObject,'Value'); 
% Update handles structure
guidata(hObject, handles);

% % --- Executes on selection change in list_selchan.
% function list_selchan_Callback(hObject, eventdata, handles)
% % hObject    handle to list_selchan (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: contents = cellstr(get(hObject,'String')) returns list_selchan contents as cell array
% %        contents{get(hObject,'Value')} returns selected item from list_selchan
% 
% idx_selchan = get(hObject,'Value');
% switch idx_selchan
%     case 1
%         handles.powchan = find(strcmp(upper(chantype(handles.Dmeg{1})),'EEG'));
%     case 2
%         [filename, pathname] = uigetfile('.txt','Open the Selection');
%         namefile = fullfile(pathname,filename);
%         fid = fopen (namefile,'r');
%         file = fread(fid,'uint8=>char')';
%         fclose(fid);
%         i = 1;
%         handles.powchan = cell(0);
%         while i<length(file)
%             init = i;
%             while ~strcmpi(file(i),' ') && i<length(file)      
%                 i=i+1;
%             end 
%             handles.powchan = [handles.powchan file(init : i-1)];
%             i=i+1;   
%         end
%     case 3
%         fprintf('to be done')
% end
% % Update handles structure
% guidata(hObject, handles);

% % --- Executes during object creation, after setting all properties.
% function list_selchan_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to list_selchan (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: listbox controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end



function sel_epoch_Callback(hObject, eventdata, handles)
% hObject    handle to sel_epoch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sel_epoch as text
%        str2double(get(hObject,'String')) returns contents of sel_epoch as a double
handles.epoch = str2double(get(hObject,'String'));
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function sel_epoch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sel_epoch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_process.
function push_process_Callback(hObject, eventdata, handles)
% hObject    handle to push_process (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% send data to the different process steps in the following order

[pathn id ext] = fileparts(handles.file);

fprintf(1,'===========================================\n');
fprintf(1,'TREATING SUBJECT %s \n',id);
fprintf(1,'===========================================\n');

if handles.flags_reref
    if handles.avgmas
        handles.Dmeg = csg_reref(handles);
        [pathn filen ext] = fileparts(handles.file);
        handles.file = fullfile(pathn,['M' id ext]);
        fprintf('Rereferenced file is: %s \n',char(['M' filen ext]));
    end
end 

if handles.flags_badchan
    handles.Dmeg = csg_badchannels(handles);
end 

if handles.flags_interpol
    handles.Dmeg = csg_interpol(handles);
end

if handles.flags_badepoch
    handles.Dmeg = CSG_artefact(handles);
end

if handles.flags_ps
    handles.plot = 1;
    handles.Dmeg = csg_powerspect(handles);
end 
delete(handles.figure1);


% --- Executes on button press in sel_psgood.
function sel_psgood_Callback(hObject, eventdata, handles)
% hObject    handle to sel_psgood (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sel_psgood
handles.flags_ps = get(hObject,'Value');
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in sel_pswholegoods.
function sel_pswholegoods_Callback(hObject, eventdata, handles)
% hObject    handle to sel_pswholegoods (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sel_pswholegoods


% --- Executes on button press in sel_psselchan.
function sel_psselchan_Callback(hObject, eventdata, handles)
% hObject    handle to sel_psselchan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sel_psselchan


% --- Executes on button press in push_interpol.
function push_interpol_Callback(hObject, eventdata, handles)
% hObject    handle to push_interpol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of push_interpol
handles.flags_interpol = get(hObject,'Value');
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in sel_badepoch.
function sel_badepoch_Callback(hObject, eventdata, handles)
% hObject    handle to sel_badepoch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sel_badepoch
handles.flags_badepoch = get(hObject,'Value');
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in reref_2mastoid.
function reref_2mastoid_Callback(hObject, eventdata, handles)
% hObject    handle to reref_2mastoid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of reref_2mastoid
handles.flags_reref = get(hObject,'Value');
% Update handles structure
guidata(hObject, handles);
