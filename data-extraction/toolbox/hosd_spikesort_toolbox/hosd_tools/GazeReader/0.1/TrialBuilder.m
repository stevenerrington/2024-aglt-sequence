function varargout = TrialBuilder(varargin)
% TRIALBUILDER M-file for TrialBuilder.fig
%      TRIALBUILDER, by itself, creates a new TRIALBUILDER or raises the existing
%      singleton*.
%
%      H = TRIALBUILDER returns the handle to a new TRIALBUILDER or the handle to
%      the existing singleton*.
%
%      TRIALBUILDER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRIALBUILDER.M with the given input arguments.
%
%      TRIALBUILDER('Property','Value',...) creates a new TRIALBUILDER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TrialBuilder_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TrialBuilder_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

% Edit the above text to modify the response to help TrialBuilder

% Last Modified by GUIDE v2.5 27-Dec-2007 13:25:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TrialBuilder_OpeningFcn, ...
                   'gui_OutputFcn',  @TrialBuilder_OutputFcn, ...
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


% --- Executes just before TrialBuilder is made visible.
function TrialBuilder_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TrialBuilder (see VARARGIN)

% Choose default command line output for TrialBuilder
handles.output = hObject;

setappdata(handles.figure1,'parent',varargin{1});
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TrialBuilder wait for user response (see UIRESUME)
% uiwait(handles.figure1);

Update(hObject, eventdata, handles)
% --- Outputs from this function are returned to the command line.
function varargout = TrialBuilder_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in eventList.
function eventList_Callback(hObject, eventdata, handles)
% hObject    handle to eventList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns eventList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from eventList


% --- Executes during object creation, after setting all properties.
function eventList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eventList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in trialList.
function trialList_Callback(hObject, eventdata, handles)
% hObject    handle to trialList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns trialList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from trialList


% --- Executes during object creation, after setting all properties.
function trialList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trialList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function offset_Callback(hObject, eventdata, handles)
% hObject    handle to offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of offset as text
%        str2double(get(hObject,'String')) returns contents of offset as a double


% --- Executes during object creation, after setting all properties.
function offset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function markTrialStart_Callback(hObject, eventdata, handles)
% hObject    handle to markTrialStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 


parent = getappdata(handles.figure1,'parent');
expEventData = getappdata(parent,'expEventData');
% trialData = getappdata(parent,'trialData');

% sttimes= [trialData.trials.startTime];
% stptimes= [trialData.stopTime];

selected = get(handles.eventList,'value');
selected(selected == 1) = [];

if length(selected) ~=1
    return
end
selected = get(handles.eventList,'value');
selxdat = expEventData.xdat(selected-1);

times = [selxdat.startT]-expEventData.xdat(1).startT;


if isfield(expEventData,'trialOnsets')
    oldto = expEventData.trialOnsets;
    oldtoc = expEventData.trialOnsetCodes;
    
    drop = ismember(times,oldto);
    selxdat(drop) = [];
    
else
    oldto = [];
    oldtoc = [];
end

expEventData.trialOnsets = cat(2,oldto,[selxdat.startT]-expEventData.xdat(1).startT)+ str2double(get(handles.offset,'string'));
expEventData.trialOnsetCodes = cat(2,oldtoc,[selxdat.id]);
    
[q,srt] = sort(expEventData.trialOnsets);
expEventData.trialOnsets =  expEventData.trialOnsets(srt);
expEventData.trialOnsetCodes =  expEventData.trialOnsetCodes(srt);


setappdata(parent,'expEventData',expEventData);
Update(hObject, eventdata, handles)


% --------------------------------------------------------------------
function markTrialStop_Callback(hObject, eventdata, handles)
% hObject    handle to markTrialStop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


parent = getappdata(handles.figure1,'parent');
expEventData = getappdata(parent,'expEventData');
% trialData = getappdata(parent,'trialData');

% sttimes= [trialData.trials.startTime];
% stptimes= [trialData.stopTime];

selected = get(handles.eventList,'value');
selected(selected == 1) = [];

if length(selected) ~=1
    return
end

selected = get(handles.eventList,'value');

selxdat = expEventData.xdat(selected-1);

times = [selxdat.startT]-expEventData.xdat(1).startT;


if isfield(expEventData,'trialEnds')
    oldte = expEventData.trialEnds;
    oldtec = expEventData.trialEndCodes;
    
    drop = ismember(times,oldte);
    selxdat(drop) = [];
    
else
    oldte = [];
    oldtec = [];
end

expEventData.trialEnds = cat(2,oldte,[selxdat.startT]-expEventData.xdat(1).startT) + str2double(get(handles.offset,'string'));
expEventData.trialEndCodes = cat(2,oldtec,[selxdat.id]);
    
[q,srt] = sort(expEventData.trialEnds);
expEventData.trialEnds =  expEventData.trialEnds(srt);
expEventData.trialEndCodes =  expEventData.trialEndCodes(srt);


setappdata(parent,'expEventData',expEventData);
Update(hObject, eventdata, handles)


% --------------------------------------------------------------------
function eventMenu_Callback(hObject, eventdata, handles)
% hObject    handle to eventMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function markCodeTrialStart_Callback(hObject, eventdata, handles)
% hObject    handle to markCodeTrialStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


parent = getappdata(handles.figure1,'parent');
expEventData = getappdata(parent,'expEventData');
% trialData = getappdata(parent,'trialData');

% sttimes= [trialData.trials.startTime];
% stptimes= [trialData.stopTime];

selected = get(handles.eventList,'value');
selected(selected == 1) = [];

if length(selected) ~=1
    return
end

codes = [expEventData.xdat.id];

selected = find( codes == expEventData.xdat(selected-1).id);

selxdat = expEventData.xdat(selected);

times = [selxdat.startT]-expEventData.xdat(1).startT;


if isfield(expEventData,'trialOnsets')
    oldto = expEventData.trialOnsets;
    oldtoc = expEventData.trialOnsetCodes;
    
    drop = ismember(times,oldto);
    selxdat(drop) = [];
    
else
    oldto = [];
    oldtoc = [];
end

expEventData.trialOnsets = cat(2,oldto,[selxdat.startT]-expEventData.xdat(1).startT)+ str2double(get(handles.offset,'string'));
expEventData.trialOnsetCodes = cat(2,oldtoc,[selxdat.id]);
    
[q,srt] = sort(expEventData.trialOnsets);
expEventData.trialOnsets =  expEventData.trialOnsets(srt);
expEventData.trialOnsetCodes =  expEventData.trialOnsetCodes(srt);


setappdata(parent,'expEventData',expEventData);
Update(hObject, eventdata, handles)

% --------------------------------------------------------------------
function markCodeTrialEnd_Callback(hObject, eventdata, handles)
% hObject    handle to markCodeTrialEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



parent = getappdata(handles.figure1,'parent');
expEventData = getappdata(parent,'expEventData');
% trialData = getappdata(parent,'trialData');

% sttimes= [trialData.trials.startTime];
% stptimes= [trialData.stopTime];

selected = get(handles.eventList,'value');
selected(selected == 1) = [];

if length(selected) ~=1
    return
end

codes = [expEventData.xdat.id];

selected = find( codes == expEventData.xdat(selected-1).id);

selxdat = expEventData.xdat(selected);

times = [selxdat.startT]-expEventData.xdat(1).startT;


if isfield(expEventData,'trialEnds')
    oldte = expEventData.trialEnds;
    oldtec = expEventData.trialEndCodes;
    
    drop = ismember(times,oldte);
    selxdat(drop) = [];
    
else
    oldte = [];
    oldtec = [];
end

expEventData.trialEnds = cat(2,oldte,[selxdat.startT]-expEventData.xdat(1).startT) + str2double(get(handles.offset,'string'));
expEventData.trialEndCodes = cat(2,oldtec,[selxdat.id]);
    
[q,srt] = sort(expEventData.trialEnds);
expEventData.trialEnds =  expEventData.trialEnds(srt);
expEventData.trialEndCodes =  expEventData.trialEndCodes(srt);


setappdata(parent,'expEventData',expEventData);
Update(hObject, eventdata, handles)


% --------------------------------------------------------------------
function Update(hObject,eventData,handles)


UpdateTrialData(hObject,eventData,handles)

parent = getappdata(handles.figure1,'parent');

expEventData = getappdata(parent,'expEventData');

currentDataSet = getappdata(parent,'CurrentDataSet');

trialData = getappdata(parent,'trialData');

if isempty(currentDataSet)
    currentDataSet = 1;
    setappdata(parent,'CurrentDataSet',currentDataSet )
end


expEventData = expEventData(currentDataSet);

if isempty(expEventData)
    error('There appears to be no loaded data set.')
end
%String for the event Data list box
evtcodenum = cellfun(@num2str,{expEventData.xdat.id},'uniformoutput',0);
evtcodes = cellfun(@num2str,expEventData.codes([expEventData.xdat.id]),'uniformoutput',0);
cstrlen = cellfun(@length,evtcodenum);
spaces1 = cellfun(@(n) repmat(' ',1,n),mat2cell(10-2*cstrlen,1,ones(1,length(cstrlen))),'UniformOutput',false);
evt = [expEventData.xdat.startT];
evt = evt-evt(1);
evttimes = cellfun(@num2str,mat2cell(evt,1,ones(1,length(evt))),'uniformoutput',0);
tstrlen = cellfun(@length,evttimes);
spaces2 = cellfun(@(n) repmat(' ',1,n),mat2cell(15-2*tstrlen,1,ones(1,length(cstrlen))),'UniformOutput',false);
liststr = strcat(evttimes,spaces2,evtcodenum,spaces1,evtcodes);

liststr = [{'Time (samp)    code#     code'},liststr];

set(handles.eventList,'string',liststr)


%String for the trial Data list box
if ~isempty([trialData.trials.number])
    
    trnum = cellfun(@num2str,{trialData.trials.number},'uniformoutput',0);
    trlen = cellfun(@length,trnum);
    spaces1 = cellfun(@(n) repmat(' ',1,n),mat2cell(8-2*trlen,1,ones(1,length(trlen))),'UniformOutput',false);

    startTime = cellfun(@num2str,{trialData.trials.startTime},'uniformoutput',0);
    sttlen = cellfun(@length,startTime);
    spaces2 = cellfun(@(n) repmat(' ',1,n),mat2cell(20-2*sttlen,1,ones(1,length(sttlen))),'UniformOutput',false);

    startcodes = cellfun(@num2str,expEventData.codes([trialData.trials.startCode]),'uniformoutput',0);
    startlen = cellfun(@length,startcodes);
    spaces3 = cellfun(@(n) repmat(' ',1,n),mat2cell(30-2*startlen,1,ones(1,length(startlen))),'UniformOutput',false);

    stopTime = cellfun(@num2str,{trialData.trials.stopTime},'uniformoutput',0);
    stoptlen = cellfun(@length,stopTime );
    spaces4 = cellfun(@(n) repmat(' ',1,n),mat2cell(20-2*stoptlen ,1,ones(1,length(stoptlen))),'UniformOutput',false);

    stopcodes = cellfun(@num2str,expEventData.codes([trialData.trials.stopCode]),'uniformoutput',0);
    stoplen = cellfun(@length,stopcodes);
    spaces5 = cellfun(@(n) repmat(' ',1,n),mat2cell(30-2*stoplen,1,ones(1,length(stoplen))),'UniformOutput',false);

    trliststr = strcat(trnum,spaces1,startTime,spaces2,startcodes,spaces3,stopTime,spaces4,stopcodes,spaces5);
else
    trliststr = {};
end
trliststr = [{'Trial#        StartT       start code        StopT        stop code'},trliststr ];

set(handles.trialList,'string',trliststr )


% fg = figure;
% set(fg,'MenuBar','none');
% lb = uicontrol(fg,'String',liststr,'Style','Listbox','units','normalized','string',liststr,'position',[.1 .1 .8 .8]);


%-------------------------------------------

function UpdateTrialData(hObject,eventData,handles)

parent = getappdata(handles.figure1,'parent');

trialData = getappdata(parent,'trialData');
currentDataSet = getappdata(parent,'CurrentDataSet');

expEventData = getappdata(parent,'expEventData');

if isempty(expEventData)
    error('No data set is loaded yet!')
end

evt = expEventData(currentDataSet); 

if ~isfield(evt,'trialEnds') || ~isfield(evt,'trialOnsets') || length(evt.trialEnds) ~= length(evt.trialOnsets) 
    return
end
  
trialdat = makeTrialData('startTime',evt.trialOnsets,'stopTime', evt.trialEnds,'startCode',...
    evt.trialOnsetCodes,'stopCode',evt.trialEndCodes,'code',(1:length(evt.trialOnsets))+ trialData(currentDataSet).codeincr);

trialData(currentDataSet).trials = trialdat; 
trialData(currentDataSet).codeincr = max([trialdat.code]);

setappdata(parent,'trialData',trialData);

tmfuns = getappdata(parent,'trialManagerFunctions');
tmhandle = getappdata(parent,'trialManager');

if ishandle(tmhandle)
    tmfuns.update;
end
