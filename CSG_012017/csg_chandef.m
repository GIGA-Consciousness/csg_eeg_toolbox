function varargout = csg_chandef(varargin)
% CSG_CHANDEF MATLAB code for csg_chandef.fig
%      CSG_CHANDEF, by itself, creates a new CSG_CHANDEF or raises the existing
%      singleton*.
%
%      H = CSG_CHANDEF returns the handle to a new CSG_CHANDEF or the handle to
%      the existing singleton*.
%
%      CSG_CHANDEF('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CSG_CHANDEF.M with the given input arguments.
%
%      CSG_CHANDEF('Property','Value',...) creates a new CSG_CHANDEF or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before csg_chandef_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to csg_chandef_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help csg_chandef

% Last Modified by GUIDE v2.5 02-Nov-2016 10:12:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @csg_chandef_OpeningFcn, ...
                   'gui_OutputFcn',  @csg_chandef_OutputFcn, ...
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


% --- Executes just before csg_chandef is made visible.
function csg_chandef_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to csg_chandef (see VARARGIN)

% Choose default command line output for csg_chandef
set(0,'CurrentFigure',handles.figure1);
handles.output = hObject;
load CRC_electrodes.mat;
handles.names     = names;
handles.pos       = pos';
handles.crc_types = crc_types;

if isempty(varargin) || ~isfield(varargin{1},'file')
    % Filter for vhdr, mat and edf files
    prefile = spm_select(1, 'any', 'Select imported EEG file','' ...
            ,pwd,'\.[mMvVeErR][dDhHaA][fFDdTtwW]');
    handles.Dmeg{1} = crc_eeg_load(deblank(prefile));
    file = fullfile(handles.Dmeg{1}.path,handles.Dmeg{1}.fname);
    handles.file{1} = file;
    handles.chan{1} = upper(chanlabels(handles.Dmeg{1}));
else
    handles.file = varargin{1}.file;
    prefile = deblank(handles.file);
    index = varargin{1}.index;
    handles.Dmeg{1} = varargin{1}.Dmeg{1};
    if isempty(index)
        index=1:nchannels(handles.Dmeg{1});
    end
    % fill the list of channels
    set(handles.list_available,'String',upper(chanlabels(handles.Dmeg{1},varargin{1}.index)));
    diff    =   setdiff(upper(chanlabels(handles.Dmeg{1})),upper(chanlabels(handles.Dmeg{1},varargin{1}.index)));
    set(handles.list_selected,'String',diff);
    % display selected channels
    [dumb1,dumb2,index2]    =   intersect(upper(chanlabels(handles.Dmeg{1},varargin{1}.index)),upper(handles.names));

    idxred=index2(find(handles.crc_types(index2)<-1));
    idxblue=index2(find(handles.crc_types(index2)>-2));

    xred=handles.pos(1,idxred);
    yred=handles.pos(2,idxred);

    xblu=handles.pos(1,idxblue);
    yblu=handles.pos(2,idxblue);

    cleargraph(handles)
    hold on
    plot(xred,yred,'r+'), plot(xblu,yblu,'b+')
    hold off
    if and(length(xblu)==0,length(xred)==0)
        cleargraph(handles)
    end

    xlim([0 1])
    ylim([0 1])
    handles.chan{1}=get(handles.Select,'String');
end

handles.eeg = [meegchannels(handles.Dmeg{1},'EEG') meegchannels(handles.Dmeg{1},'LFP')];
handles.meg = meegchannels(handles.Dmeg{1},'MEG');
handles.eog = eogchannels(handles.Dmeg{1});
handles.emg = emgchannels(handles.Dmeg{1});
handles.ecg = ecgchannels(handles.Dmeg{1});
handles.other = setdiff(1:size(handles.Dmeg{1},1),[handles.eeg, handles.meg, handles.eog, handles.emg, handles.ecg]);
handles.typechan = chantype(handles.Dmeg{1});
handles.labelchan = chanlabels(handles.Dmeg{1});

poptype = {'EEG','MEG','EOG','EMG','ECG','REF','Other'};
chanset = handles.chan{1};
% display list of channels
set(handles.list_available,'String',chanset);
set(handles.chantype,'String',poptype,'Value',length(poptype));

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes csg_dis_selchan wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% UIWAIT makes csg_chandef wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = csg_chandef_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in list_available.
function list_available_Callback(hObject, eventdata, handles)
% hObject    handle to list_available (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns list_available contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_available
contents = get(hObject,'String');
if length(contents)==0
else
    %Remove the "activated" item from the list "Available Channels"
    [dumb1,dumb2,index]=intersect(contents{get(hObject,'Value')},contents);
    temp=[contents(1:index-1) ; contents(index+1:length(contents))];
    set(handles.list_available,'String',temp);
     
    %Add the "activated" in the list "Selected Channels"
    if length(get(handles.list_selected,'String'))==0
        temp={contents{get(hObject,'Value')}};
    else
        temp=[contents{get(hObject,'Value')} ; get(handles.list_selected,'String')];
    end
    set(handles.list_selected,'String',temp);

    %Prevent crashing if the first/last item of the list is selected.
    set(handles.list_available,'Value',max(index-1,1));
    set(handles.list_selected,'Value',1);

    [dumb1,dumb2,index]=intersect(upper(temp),upper(handles.names));

    idxred=index(find(handles.crc_types(index)<-1));
    idxblue=index(find(handles.crc_types(index)>-2));

    xred=handles.pos(1,idxred);
    yred=handles.pos(2,idxred);

    xblu=handles.pos(1,idxblue);
    yblu=handles.pos(2,idxblue);

    cleargraph(handles)
    hold on
    plot(xred,yred,'r+'), plot(xblu,yblu,'b+')
    hold off
    if and(length(xblu)==0,length(xred)==0)
        cleargraph(handles)
    end

    xlim([0 1])
    ylim([0 1])

    set(handles.axes1,'XTick',[]);
    set(handles.axes1,'YTick',[]);

    %if multiple comparison, no more than one channel
    contents=get(handles.list_selected,'String');
    if isfield(handles,'multcomp') && (handles.multcomp && length(contents)>1 || isempty(contents))
        set(handles.PLOT,'enable','off');
        set(handles.PLOT,'ForegroundColor',[1 0 0]);
        beep
        disp('Select only one channel for multiple files comparison')
    elseif isfield(handles,'multcomp') && (handles.multcomp && length(contents)==1)
        set(handles.PLOT,'enable','on');
        set(handles.PLOT,'ForegroundColor',[0 0 0]);
    end
end
% Indicate this selection doesn't come from an selection uploaded
set(handles.push_save,'value',1);
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function list_available_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_available (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in list_selected.
function list_selected_Callback(hObject, eventdata, handles)
% hObject    handle to list_selected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns list_selected contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_selected
contents = get(hObject,'String');
if isempty(contents)
else
    %Remove the "activated" item from the list "Available Channels"
    [dumb1,dumb2,index] = intersect(contents{get(hObject,'Value')},contents);
    temp = [contents(1:index-1) ; contents(index+1:length(contents))];
    set(handles.list_selected,'String',temp);

    [dumb1,dumb2,index2]=intersect(upper(temp),upper(handles.names));

    idxred=index2(find(handles.crc_types(index2)<-1));
    idxblue=index2(find(handles.crc_types(index2)>-2));

    xred=handles.pos(1,idxred);
    yred=handles.pos(2,idxred);

    xblu=handles.pos(1,idxblue);
    yblu=handles.pos(2,idxblue);

    cleargraph(handles)
    hold on
    plot(xred,yred,'r+'), plot(xblu,yblu,'b+')
    hold off
    if and(length(xblu)==0,length(xred)==0)
        cleargraph(handles)
    end

    xlim([0 1])
    ylim([0 1])

    set(handles.axes1,'XTick',[]);
    set(handles.axes1,'YTick',[]);

    %Add the "activated" in the list "Selected Channels"
    if length(get(handles.list_available,'String'))==0
        temp={contents{get(hObject,'Value')}};
    else
        temp=[contents{get(hObject,'Value')} ; get(handles.list_available,'String')];
    end
    set(handles.list_available,'String',temp);

    %Prevent crashing if the first/last item of the list is selected.
    set(handles.list_selected,'Value',max(index-1,1));
    set(handles.list_available,'Value',1);

    %if multiple comparison, no more than one channel
    contents=get(handles.list_selected,'String');
    if isfield(handles,'multcomp') && (handles.multcomp && length(contents)>1 || isempty(contents))
        set(handles.PLOT,'enable','off');
        set(handles.PLOT,'ForegroundColor',[1 0 0]);
        beep
        disp('Select only one channel for multiple files comparison')
    elseif isfield(handles,'multcomp') && (handles.multcomp && length(contents)==1)
        set(handles.PLOT,'enable','on');
        set(handles.PLOT,'ForegroundColor',[0 0 0]);
    end
end
% Indicate this selection doesn't come from an selection uploaded
set(handles.push_save,'value',1);
% set(handles.load,'Value',0);
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function list_selected_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_selected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_save.
function push_save_Callback(hObject, eventdata, handles)
% hObject    handle to push_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% update the chantype 
typeselect = get(handles.chantype,'String');
chanselect = get(handles.list_selected,'String');
[gum idxchan] = intersect(chanlabels(handles.Dmeg{1}),chanselect);
handles.typechan(idxchan) = typeselect(get(handles.chantype,'Value'));
handles.Dmeg{1} = chantype(handles.Dmeg{1},1:size(handles.Dmeg{1},1),handles.typechan);
% Indicate this selection is saved
save(handles.Dmeg{1});
set(handles.push_save,'value',0);
guidata(hObject, handles);

% --- Executes on button press in push_exit.
function push_exit_Callback(hObject, eventdata, handles)
% hObject    handle to push_exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
varargout{1} = handles.output;
delete(handles.figure1)

% --- Executes on selection change in chantype.
function chantype_Callback(hObject, eventdata, handles)
% hObject    handle to chantype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns chantype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from chantype
contents = get(hObject,'String');
if isempty(contents)
else
    typeselected = contents(get(hObject,'Value'));    
    idxchan = find(strcmp(chantype(handles.Dmeg{1}),typeselected));
    idx_available = setdiff(1:size(handles.Dmeg{1},1),idxchan);
    set(handles.list_selected,'String',chanlabels(handles.Dmeg{1},idxchan));
    set(handles.list_available,'String',chanlabels(handles.Dmeg{1},idx_available));
    
    [dumb1,dumb2,index2]=intersect(upper(chanlabels(handles.Dmeg{1},idxchan)),upper(handles.names));

    idxred=index2(find(handles.crc_types(index2)<-1));
    idxblue=index2(find(handles.crc_types(index2)>-2));

    xred=handles.pos(1,idxred);
    yred=handles.pos(2,idxred);

    xblu=handles.pos(1,idxblue);
    yblu=handles.pos(2,idxblue);

    cleargraph(handles)
    hold on
    plot(xred,yred,'r+'), plot(xblu,yblu,'b+')
    hold off
    if and(length(xblu)==0,length(xred)==0)
        cleargraph(handles)
    end

    xlim([0 1])
    ylim([0 1])

    set(handles.axes1,'XTick',[]);
    set(handles.axes1,'YTick',[]);
end 
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function chantype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chantype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in select_all.
function select_all_Callback(hObject, eventdata, handles)
% hObject    handle to select_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.list_selected,'String',upper(chanlabels(handles.Dmeg{1})));
set(handles.list_available,'String',cell(0));

[dumb1,dumb2,index]=intersect(upper(chanlabels(handles.Dmeg{1})),upper(handles.names));

idxred=index(find(handles.crc_types(index)<-1));
idxblue=index(handles.crc_types(index)>-2);

xred=handles.pos(1,idxred);
yred=handles.pos(2,idxred);

xblu=handles.pos(1,idxblue);
yblu=handles.pos(2,idxblue);

cleargraph(handles)
hold on
plot(xred,yred,'r+'), plot(xblu,yblu,'b+')
hold off
if and(length(xblu)==0,length(xred)==0)
    cleargraph(handles)
end

xlim([0 1])
ylim([0 1])

set(handles.axes1,'XTick',[]);
set(handles.axes1,'YTick',[]);

% Indicate this selection doesn't come from an selection uploaded
set(handles.push_save,'value',1);
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in deselect_all.
function deselect_all_Callback(hObject, eventdata, handles)
% hObject    handle to deselect_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.list_available,'String',upper(chanlabels(handles.Dmeg{1})));
set(handles.list_selected,'String',cell(0));

cleargraph(handles)
xlim([0 1])
ylim([0 1])
set(handles.axes1,'XTick',[]);
set(handles.axes1,'YTick',[]);
% Indicate this selection doesn't come from an selection uploaded
set(handles.push_save,'visible','on');
% set(handles.load,'Value',0);
% Update handles structure
guidata(hObject, handles);
%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
%%% subfunction %%%
%%%%%%%%%%%%%%%%%%%
function cleargraph(handles)

A=get(handles.figure1,'Children');
idx=find(strcmp(get(A,'Type'),'axes')==1);
try
    delete(get(A(idx),'Children'))
end


