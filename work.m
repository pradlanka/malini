function varargout = work(varargin)
%global handles.lmn;
% WORK MATLAB code for work.fig
%      WORK, by itself, creates a new WORK or raises the existing
%      singleton*.
%
%      H = WORK returns the handle to a new WORK or the handle to
%      the existing singleton*.
%
%      WORK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WORK.M with the given input arguments.
%
%      WORK('Property','Value',...) creates a new WORK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before work_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to work_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help work

% Last Modified by GUIDE v2.5 26-Jun-2020 17:57:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @work_OpeningFcn, ...
                   'gui_OutputFcn',  @work_OutputFcn, ...
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


% --- Executes just before work is made visible.
function work_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to work (see VARARGIN)

% Choose default command line output for work
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
% UIWAIT makes work wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = work_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Filepath=pwd;
folder_name=get(hObject,'string');
handles.folder_name=folder_name;
guidata(hObject,handles);
assignin('base','folder_name',handles.folder_name);



% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse.
function browse_Callback(hObject, eventdata, handles)
% hObject    handle to browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Filepath=pwd;
data_folder_name=uigetdir('c:\');
handles.data_folder_name=data_folder_name;
set(handles.edit1,'String',data_folder_name);
%cd(data_folder_name);
%Groupdir=dir;
%for movin=1:length(Groupdir)-2
%mkdir(fullfile(folder_name,'Preprocessed-data',Groupdir(movin+2).name));
%movefile(Groupdir(movin+2).name,'Preprocessed-data');
%end
cd(Filepath);
%handles.folder_name = folder_name;
%guidata(hObject,handles);
%addpath(genpath(folder_name));
assignin('base','data_folder_name',handles.data_folder_name);


%[filename pathname] = uigetfile({'.'}, 'File selector');
%fullpathname = strcat(pathname,filename);
%text=fileread(fullpathname);
%set(handles.edit1,'string',fullpathname);
%assignin('base','folder',[filename pathname])




function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in adhd.
function adhd_Callback(hObject, eventdata, handles)
% hObject    handle to adhd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in controls.
function controls_Callback(hObject, eventdata, handles)
% hObject    handle to controls (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in controls1.
function controls1_Callback(hObject, eventdata, handles)
% hObject    handle to controls1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in staticfc.
function staticfc_Callback(hObject, eventdata, handles)
% hObject    handle to staticfc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%addpath(genpath('C:\Users\xyz\Documents\MATLAB\gui'))
%a = get(handles.sfc,'Value')
%if(a == 1)
 %  Controls_together
 %PTSD_together
 %staticFC_individual
 %data_create_sfc
%end
% Hint: get(hObject,'Value') returns toggle state of staticfc


% --- Executes on button press in sec.
function sec_Callback(hObject, eventdata, handles)
% hObject    handle to sec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%addpath(genpath('C:\Users\xyz\Documents\MATLAB\gui'))
%stec = get(handles.time_series,'Value')
%if(stec == 1)
%Controls_together
 %PTSD_together
  % secout
   %data_create_sec
   % assignin('base','con',conmat)
%end
% Hint: get(hObject,'Value') returns toggle state of sec


% --- Executes on button press in dec.
function dec_Callback(hObject, eventdata, handles)
% hObject    handle to dec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dec


% --- Executes on button press in dfc.
function dfc_Callback(hObject, eventdata, handles)
% hObject    handle to dfc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%addpath(genpath('C:\Users\xyz\Documents\MATLAB\gui'))
%dyfc = get(handles.dfc,'Value')
%if(dyfc == 1)
 %Controls_together
 %PTSD_together
 %dfcout
 %data_create
  % dfcout
%end
% Hint: get(hObject,'Value') returns toggle state of dfc


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in deconvolution.
function deconvolution_Callback(hObject, eventdata, handles)
% hObject    handle to deconvolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%deco=get(handles.deconvolution,'Value')
%if(deco==1)
 %   mat_extraction
  %  perform_dec
   % after_perform_dec
    %post_decon_extract_timeseries
%end
% Hint: get(hObject,'Value') returns toggle state of deconvolution


% --- Executes on button press in recursivecluster.
function recursivecluster_Callback(hObject, eventdata, handles)
% hObject    handle to recursivecluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of recursivecluster


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function initialcluster_Callback(hObject, eventdata, handles)
% hObject    handle to initialcluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%incl = str2double(get( hObject, 'String' ));
%handles.incl=incl;
%guidata(hObject,handles);
%assignin('base','incl',handles.incl);
global incl;
incl = str2num(get(hObject,'String'));
handles.incl=incl;
guidata(hObject, handles);
assignin('base','incl',handles.incl);
% Hints: get(hObject,'String') returns contents of initialcluster as text
%        str2double(get(hObject,'String')) returns contents of initialcluster as a double


% --- Executes during object creation, after setting all properties.
function initialcluster_CreateFcn(hObject, eventdata, handles)
% hObject    handle to initialcluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function finalcluster_Callback(hObject, eventdata, handles)
% hObject    handle to finalcluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%fincl = str2double(get( hObject, 'String' ));
%handles.fincl=fincl;
%guidata(hObject,handles);
%assignin('base','fincl',handles.fincl);
global fincl;
fincl = str2num(get(hObject,'String'));
handles.fincl=fincl;
guidata(hObject, handles);
assignin('base','fincl',handles.fincl);
% Hints: get(hObject,'String') returns contents of finalcluster as text
%        str2double(get(hObject,'String')) returns contents of finalcluster as a double


% --- Executes during object creation, after setting all properties.
function finalcluster_CreateFcn(hObject, eventdata, handles)
% hObject    handle to finalcluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function decreasecluster_Callback(hObject, eventdata, handles)
% hObject    handle to decreasecluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global deccl;
deccl = str2double(get( hObject, 'String' ));
handles.deccl=deccl;
guidata(hObject,handles);
assignin('base','deccl',handles.deccl);
% Hints: get(hObject,'String') returns contents of decreasecluster as text
%        str2double(get(hObject,'String')) returns contents of decreasecluster as a double


% --- Executes during object creation, after setting all properties.
function decreasecluster_CreateFcn(hObject, eventdata, handles)
% hObject    handle to decreasecluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function endpath_Callback(hObject, eventdata, handles)
% hObject    handle to endpath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global epat;
epat = str2double(get( hObject, 'String' ));
handles.epat=epat;
guidata(hObject,handles);
assignin('base','epat',handles.epat);
% Hints: get(hObject,'String') returns contents of endpath as text
%        str2double(get(hObject,'String')) returns contents of endpath as a double


% --- Executes during object creation, after setting all properties.
function endpath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to endpath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function folds_Callback(hObject, eventdata, handles)
% hObject    handle to folds (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
nofo = str2double(get( hObject, 'String' ));
handles.nofo=nofo;
guidata(hObject,handles);
assignin('base','nofo',handles.nofo);
% Hints: get(hObject,'String') returns contents of folds as text
%        str2double(get(hObject,'String')) returns contents of folds as a double


% --- Executes during object creation, after setting all properties.
function folds_CreateFcn(hObject, eventdata, handles)
% hObject    handle to folds (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function resampling_Callback(hObject, eventdata, handles)
% hObject    handle to resampling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
resa = str2double(get( hObject, 'String' ))
handles.resa=resa;
guidata(hObject,handles);
assignin('base','resa',handles.resa);
% Hints: get(hObject,'String') returns contents of resampling as text
%        str2double(get(hObject,'String')) returns contents of resampling as a double


% --- Executes during object creation, after setting all properties.
function resampling_CreateFcn(hObject, eventdata, handles)
% hObject    handle to resampling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function classes_Callback(hObject, eventdata, handles)
% hObject    handle to classes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of classes as text
%        str2double(get(hObject,'String')) returns contents of classes as a double
nocl = str2double(get( hObject, 'String' ))
handles.nocl=nocl;
guidata(hObject,handles);
assignin('base','nocl',handles.nocl);

% --- Executes during object creation, after setting all properties.
function classes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to classes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in classifier.
function classifier_Callback(hObject, eventdata, handles)
% hObject    handle to classifier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns classifier contents as cell array
%        contents{get(hObject,'Value')} returns selected item from classifier


% --- Executes during object creation, after setting all properties.
function classifier_CreateFcn(hObject, eventdata, handles)
% hObject    handle to classifier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in time_series.
function time_series_Callback(hObject, eventdata, handles)
% hObject    handle to time_series (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
presdir=pwd;
tmse=get(handles.time_series,'Value');
if(tmse==1)
    templa=strcat(presdir,'\Templates\cc200_91x109x91.nii');
    handles.templa=templa;
guidata(hObject,handles);
assignin('base','templa',handles.templa);
end
    if(tmse==2)
         templa=strcat(presdir,'\Templates\CC400.nii');
    handles.templa=templa;
guidata(hObject,handles);
assignin('base','templa',handles.templa);
    end
        if(tmse==3)
            templa=strcat(presdir,'\Templates\dos160_roi_atlas.nii');
    handles.templa=templa;
guidata(hObject,handles);
assignin('base','templa',handles.templa);
        end
if(tmse==4)
   % templa=uigetfile({'.nii'});
    %handles.templa=templa;
    [temfilnme, tempthnme]=uigetfile({'.nii'});
templa=strcat(tempthnme,temfilnme);
handles.templa=templa;
guidata(hObject,handles);
assignin('base','templa',handles.templa);
end
% Hints: contents = cellstr(get(hObject,'String')) returns time_series contents as cell array
%        contents{get(hObject,'Value')} returns selected item from time_series
%addpath(genpath('C:\Users\xyz\Documents\MATLAB\gui'))

%a = get(handles.time_series,'Value')
%if(a == 1)
 %   msgbox('wait for couple of minutes')
  % Extract_timeseries
   %msgbox('extraction done, choose sfc, dec, sec, dec')
   
%end

% --- Executes during object creation, after setting all properties.
function time_series_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_series (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in sfc.
function sfc_Callback(hObject, eventdata, handles)
% hObject    handle to sfc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%addpath(genpath('C:\Users\xyz\Documents\MATLAB\gui'))
%stfc = get(handles.sfc,'Value')
%if(stfc == 1)
%Controls_together
 %PTSD_together
  %staticFC_individual
  %data_create_sfc
   %assignin('base','result',staticFC1)
%end
% Hint: get(hObject,'Value') returns toggle state of sfc



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in dsfc.
function dsfc_Callback(hObject, eventdata, handles)
% hObject    handle to dsfc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%deconsfc=get(handles.dsfc,'Value');
%if(deconsfc==1)
 %   Controls_together_decon
  %  PTSD_together_decon
   % staticFC_decon_individual
    %data_create_decon_sfc
    
%end
% Hint: get(hObject,'Value') returns toggle state of dsfc


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5


% --- Executes on button press in ddfc.
function ddfc_Callback(hObject, eventdata, handles)
% hObject    handle to ddfc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%decondfc=get(handles.ddfc,'Value')
%if(decondfc==1)
 %   Controls_together_decon
  %  PTSD_together_decon
   % dfcout_decon
    %data_create_decon_dfc
    
    
%end
    
% Hint: get(hObject,'Value') returns toggle state of ddfc


% --- Executes on button press in dsec.
function dsec_Callback(hObject, eventdata, handles)
% hObject    handle to dsec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%deconsec=get(handles.dsec,'Value')
%if(deconsec==1)
 %   Controls_together_decon
  %  PTSD_together_decon
   % secout_decon
    %data_create_sec_decon
%end
% Hint: get(hObject,'Value') returns toggle state of dsec


% --- Executes on button press in SVMradio.
function SVMradio_Callback(hObject, eventdata, handles)
% hObject    handle to SVMradio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SVMradio


% --- Executes on button press in RCEradio.
function RCEradio_Callback(hObject, eventdata, handles)
% hObject    handle to RCEradio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RCEradio


% --- Executes when selected object is changed in rcesvm.
function rcesvm_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in rcesvm 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
blah=get(handles.rcesvm,'SelectedObject');
selecti=get(blah,'String');
switch selecti
    case 'Without RCE'
        set(handles.uibuttongroup5,'Visible','Off');
        %set(handles.checkbox13,'Enable','On');
    case 'With RCE'
        set(handles.uibuttongroup5,'Visible','On')
       % set(handles.checkbox13,'Enable','Off');
end


% --- Executes on button press in checkbox8.
function checkbox8_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.HPS=1;
dsf=get(handles.checkbox8,'Value');
if(dsf==1)

       handles.DSF = true;
else
  handles.DSF= false;
  set(handles.checkbox10,'Enable','Off');
end
guidata(hObject, handles);
assignin('base','DSF',handles.DSF);
%if((handles.HPS==true) && (handles.DSF==true))
 %   set(handles.checkbox10,'Enable','On')
%end


    
%Hint: get(hObject,'Value') returns toggle state of checkbox8


% --- Executes on button press in checkbox9.
function checkbox9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.DSF=1;
hps=get(handles.checkbox9,'Value');
if(hps==1)
    handles.HPS=true;
else
    handles.HPS=false;
    set(handles.checkbox10,'Enable','Off');
end
guidata(hObject, handles);
assignin('base','HPS',handles.HPS);

if((handles.HPS==true) && (handles.DSF==true))
    set(handles.checkbox10,'Enable','On');
elseif((handles.DSF==false) | (handles.HPS==false))
    set(handles.checkbox10,'Enable','Off');
elseif((handles.HPS==true) && (handles.DSF==true))
   set(handles.checkbox10,'Enable','On');
end
%if(handles.DSF==false)
  %   set(handles.checkbox10,'Enable','Off')
%end
%if(handles.HPS==false)
 %set(handles.checkbox10,'Enable','Off')
%end
% Hint: get(hObject,'Value') returns toggle state of checkbox9


% --- Executes on button press in checkbox10.
function checkbox10_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

madc=get(handles.checkbox10,'Value');
if(madc==1)
    handles.MADC=true;
else
    handles.MADC=false;
end
guidata(hObject, handles);
assignin('base','MADC',handles.MADC);
% Hint: get(hObject,'Value') returns toggle state of checkbox10


% --- Executes on button press in checkbox11.
function checkbox11_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%ssp=get(handles.checkbox11,'Value')
%if(ssp==1)
 %  handles.SSP=true
%else
 % handles.SSP=false
%end
%guidata(hObject, handles);
%assignin('base','SSP',handles.SSP)
% Hint: get(hObject,'Value') returns toggle state of checkbox11


% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
addpath(genpath(pwd));
Filepath=pwd;
warning('off','all');
dsf=get(handles.checkbox8,'Value');
if(dsf==1)

       handles.DSF = true;
else
     handles.DSF= false;
end
guidata(hObject, handles);
assignin('base','DSF',handles.DSF);

hps=get(handles.checkbox9,'Value');
if(hps==1)
    handles.HPS=true;
else
    handles.HPS=false;
end
guidata(hObject, handles);
assignin('base','HPS',handles.HPS);

%dsf=get(handles.checkbox8,'Value')
%hps=get(handles.checkbox9,'Value')
%if(dsf == 1 && hps==1)
 %   set(handles.checkbox10, 'Enable','On')
%end

madc=get(handles.checkbox10,'Value');
if(madc==1)
    handles.MADC=true;
else
    handles.MADC=false;
end
guidata(hObject, handles);
assignin('base','MADC',handles.MADC);

ssp=get(handles.checkbox11,'Value');
if(ssp==1)
   handles.SSP=true;
else
    handles.SSP=false;
end
guidata(hObject, handles);
assignin('base','SSP',handles.SSP);

if exist('nofo')
nofo=handles.nofo;
assignin('base','nofo',handles.nofo);
end
if exist('nocl')
nocl=handles.nocl;
assignin('base','nocl',handles.nocl);
end
if exist('resa')
resa=handles.resa;
assignin('base','resa',handles.resa);
end
%hypeelm=handles.hypeelm;
DSF=handles.DSF;
HPS=handles.HPS;
MADC=handles.MADC;
SSP=handles.SSP;
%IIP=handles.IIP;
%lmn=handles.lmn;
if exist('YourFile')
YourFile=handles.YourFile;
end
% global folder_name;
% folder_name=handles.folder_name;
global data_folder_name;
data_folder_name=get(handles.edit1,'string'); %folder_name = handles.folder_name;
global folder_name;
folder_name=data_folder_name;
%if(handles.DSF==1 && handles.HPS==1)
 %set(handles.MADS, 'Enable','On')
%end
if exist('lowv')
lowv=handles.lowv;
assignin('base','lowv',handles.lowv);
end
if exist('incred')
incred=handles.incred;
assignin('base','incred',handles.incred);
end
if exist('uppe')
uppe=handles.uppe;
assignin('base','uppe',handles.uppe);
end
global fincl incl deccl epat;

if exist('lowto')
lowto=handles.lowto;
end
if exist('incredto')
incredto=handles.incredto;
end
if exist('uppeto')
uppeto=handles.uppeto;
end



deco=get(handles.deconvolution,'Value');
 % if(deco==1|tmse==1)
  %  
%end

tmse = get(handles.time_series,'Value');
if(tmse == 1|tmse==2|tmse==3|tmse==4)
    templa=handles.templa;
    % msgbox('wait for couple of minutes')
       Extract_timeseries
   %msgbox('extraction done, choose sfc, dec, sec, dec')
end


stfc = get(handles.sfc,'Value');
if(stfc == 1)
subjects_together;
 staticFC_individual;
  data_create_sfc;
   randasp_sfc;
  global idod;
  idod='Results_sfc_';
  cool;
  clear idod;
   end
dyfc = get(handles.dfc,'Value');
if(dyfc == 1)
 subjects_together;
 dfcout;
 data_create_dfc;
 randasp_dfc;
 global idod;
   idod='Results_dfc_';
   cool;
 clearvars idod;
 end
stec=get(handles.sec,'Value');
if(stec == 1)
subjects_together;
   secout;
   data_create_sec;
   randasp_sec;
   global idod;
    idod='Results_sec_';
   cool;
   clearvars idod;
 end
  deco=get(handles.deconvolution,'Value');
if(deco==1)
     mat_extraction;
    perform_dec;
    after_perform_dec;
    post_decon_extract_timeseries;
end
   deconsec=get(handles.dsec,'Value');
  if(deconsec==1)
    subjects_together_decon;
    secout_decon;
    data_create_sec_decon;
    randasp_decon_sec;
    global idod;
  idod='Results_decon_sec_';
    cool;
    clearvars idod;
      end
deconsfc=get(handles.dsfc,'Value');
  if(deconsfc==1)
    subjects_together_decon;
    staticFC_decon_individual;
    data_create_decon_sfc;
    randasp_decon_sfc;
    global idod;
  idod='Results_decon_sfc_';
    cool;
    clearvars idod;
  end
  
 decondfc=get(handles.ddfc,'Value');
if(decondfc==1)
    subjects_together_decon;
    dfcout_decon;
    data_create_decon_dfc;
    randasp_decon_dfc;
     global idod;
  idod='Results_decon_dfc_';
    cool;
    clearvars idod;
end

loaddata = get(handles.checkbox34,'Value');
if(loaddata==1)
cool;
end
clas21 = get(handles.consencheck,'Value');
 if(clas21==1)
          global discos;
     discos='Combined_Results';
       CombinedResults_new;
        clearvars discos;
    end

   


% --- Executes during object creation, after setting all properties.
function rcesvm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rcesvm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function uibuttongroup6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uibuttongroup6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

%dsf=get(handles.checkbox8,'Value')
%hps=get(handles.checkbox9,'Value')
%if(dsf == 1 && hps==1)
 %   set(handles.checkbox10, 'Enable','on')
%end



function lowerLLimit_Callback(hObject, eventdata, handles)
% hObject    handle to lowerLLimit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%lowv=zeros(1,7);


%f(cl Hints: get(hObject,'String') returns contents of lowerLLimit as text
%        str2double(get(hObject,'String')) returns contents of lowerLLimit as a double


% --- Executes during object creation, after setting all properties.
function lowerLLimit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lowerLLimit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function increMent_Callback(hObject, eventdata, handles)
% hObject    handle to increMent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
incred = str2double(get( hObject, 'String' ));
handles.incred=incred;
guidata(hObject,handles);
assignin('base','incred',handles.incred);
% Hints: get(hObject,'String') returns contents of increMent as text
%        str2double(get(hObject,'String')) returns contents of increMent as a double


% --- Executes during object creation, after setting all properties.
function increMent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to increMent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function upperLLimit_Callback(hObject, eventdata, handles)
% hObject    handle to upperLLimit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uppe = str2double(get( hObject, 'String' ));
handles.uppe=uppe;
guidata(hObject,handles);
assignin('base','uppe',handles.uppe);
% Hints: get(hObject,'String') returns contents of upperLLimit as text
%        str2double(get(hObject,'String')) returns contents of upperLLimit as a double


% --- Executes during object creation, after setting all properties.
function upperLLimit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to upperLLimit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox13.


function secondllimit_Callback(hObject, eventdata, handles)
% hObject    handle to secondllimit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lowto = str2double(get( hObject, 'String' ));
handles.lowto=lowto;
guidata(hObject,handles);
assignin('base','lowto',handles.lowto);
% Hints: get(hObject,'String') returns contents of secondllimit as text
%        str2double(get(hObject,'String')) returns contents of secondllimit as a double


% --- Executes during object creation, after setting all properties.
function secondllimit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to secondllimit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function seconduppe_Callback(hObject, eventdata, handles)
% hObject    handle to seconduppe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
incredto = str2double(get( hObject, 'String' ));
handles.incredto=incredto;
guidata(hObject,handles);
assignin('base','incredto',handles.incredto);
% Hints: get(hObject,'String') returns contents of seconduppe as text
%        str2double(get(hObject,'String')) returns contents of seconduppe as a double


% --- Executes during object creation, after setting all properties.
function seconduppe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to seconduppe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function secondincred_Callback(hObject, eventdata, handles)
% hObject    handle to secondincred (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uppeto = str2double(get( hObject, 'String' ));
handles.uppeto=uppeto;
guidata(hObject,handles);
assignin('base','uppeto',handles.uppeto);
% Hints: get(hObject,'String') returns contents of secondincred as text
%        str2double(get(hObject,'String')) returns contents of secondincred as a double


% --- Executes during object creation, after setting all properties.
function secondincred_CreateFcn(hObject, eventdata, handles)
% hObject    handle to secondincred (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loda.
function loda_Callback(hObject, eventdata, handles)
% hObject    handle to loda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filnme, pthnme]=uigetfile({'.xlsx'});
test_path=strcat(pthnme,filnme);
handles.test_path=test_path;
%[~,~,YourFile]=xlsread(test_Path);
%handles.YourFile=YourFile;
guidata(hObject,handles);
assignin('base','test_path',handles.test_path);
%assignin('base','YourFile',handles.YourFile);
%global testdat;


% --- Executes on button press in rceknncheck.
function rceknncheck_Callback(hObject, eventdata, handles)
% hObject    handle to rceknncheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clas3 = get(handles.rceknncheck,'Value'); 
if(clas3==1)
    
     hypeknn = inputdlg({'Lower Limit','Increment By','Upper Limit','Exponential Search'},...
              'KNN -K (No of nearest neighbors)', [1 75; 1 75; 1 75; 1 75]);
handles.hypeknn=hypeknn;
guidata(hObject,handles);
assignin('base','hypeknn',handles.hypeknn);
end
% Hint: get(hObject,'Value') returns toggle state of rceknncheck


% --- Executes on button press in ldacheck.
function ldacheck_Callback(hObject, eventdata, handles)
% hObject    handle to ldacheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clas4 = get(handles.ldacheck,'Value');
if(clas4==1)
     
end
% Hint: get(hObject,'Value') returns toggle state of ldacheck


% --- Executes on button press in elmcheck.
function elmcheck_Callback(hObject, eventdata, handles)
% hObject    handle to elmcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clas2 = get(handles.elmcheck,'Value');
if(clas2==1)
      hypeelm = inputdlg({'Lower Limit','Increment By','Upper Limit','Exponential Search'},...
              'ELM- C (Regularization coefficient)', [1 75; 1 75; 1 75; 1 75]);
handles.hypeelm=hypeelm;
guidata(hObject,handles);
assignin('base','hypeelm',handles.hypeelm);
hypeelm2 = inputdlg({'Lower Limit','Increment By','Upper Limit','Exponential Search'},...
              'ELM- Gamma ', [1 75; 1 75; 1 75; 1 75]);
handles.hypeelm2=hypeelm2;
guidata(hObject,handles);
assignin('base','hypeelm2',handles.hypeelm2);
end
  % Hint: get(hObject,'Value') returns toggle state of elmcheck


% --- Executes on button press in baggedcheck.
function baggedcheck_Callback(hObject, eventdata, handles)
% hObject    handle to baggedcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clas9 = get(handles.baggedcheck,'Value');
if(clas9==1)
          
end
% Hint: get(hObject,'Value') returns toggle state of baggedcheck


% --- Executes on button press in svmcheck.
function svmcheck_Callback(hObject, eventdata, handles)
% hObject    handle to svmcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clas5 = get(handles.svmcheck,'Value');
if(clas5==1)
           hypelsvm = inputdlg({'Lower Limit','Increment By','Upper Limit','Exponential Search'},...
              'Linear-SVM- C (Regularization coefficient)', [1 75; 1 75; 1 75; 1 75]);
handles.hypelsvm=hypelsvm;
guidata(hObject,handles);
assignin('base','hypelsvm',handles.hypelsvm);
end
% Hint: get(hObject,'Value') returns toggle state of svmcheck


% --- Executes on button press in rotationcheck.
function rotationcheck_Callback(hObject, eventdata, handles)
% hObject    handle to rotationcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 clas18 = get(handles.rotationcheck,'Value');
if(clas18==1)
          
      hyperotfo = inputdlg({'Lower Limit','Increment By','Upper Limit','Exponential Search'},...
              'Rotation Forest- K (No of nearest neighbors)', [1 75; 1 75; 1 75; 1 75]);
handles.hyperotfo=hyperotfo;
guidata(hObject,handles);
assignin('base','hyperotfo',handles.hyperotfo);
 end
% Hint: get(hObject,'Value') returns toggle state of rotationcheck


% --- Executes on button press in randomcheck.
function randomcheck_Callback(hObject, eventdata, handles)
% hObject    handle to randomcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 clas16 = get(handles.randomcheck,'Value');
if(clas16==1)
            
 end
% Hint: get(hObject,'Value') returns toggle state of randomcheck


% --- Executes on button press in lvqnetcheck.
function lvqnetcheck_Callback(hObject, eventdata, handles)
% hObject    handle to lvqnetcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clas14 = get(handles.lvqnetcheck,'Value');
if(clas14==1)
             
       hypelvq = inputdlg({'Lower Limit','Increment By','Upper Limit','Exponential Search'},...
              'LVQNET- Size of hidden layer ', [1 75; 1 75; 1 75; 1 75]);
handles.hypelvq=hypelvq;
guidata(hObject,handles);
assignin('base','hypelvq',handles.hypelvq);
  end
% Hint: get(hObject,'Value') returns toggle state of lvqnetcheck


% --- Executes on button press in oknncheck.
function oknncheck_Callback(hObject, eventdata, handles)
% hObject    handle to oknncheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 %clas13 = get(handles.oknncheck,'Value');
%if(clas13==1)
  %         set(handles.uibuttongroup8,'Visible','Off');
 %     set(handles.uibuttongroup7,'Visible','On');
 % end
% Hint: get(hObject,'Value') returns toggle state of oknncheck


% --- Executes on button press in naivecheck.
function naivecheck_Callback(hObject, eventdata, handles)
% hObject    handle to naivecheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 clas6 = get(handles.naivecheck,'Value');
if(clas6==1)
            
  end
% Hint: get(hObject,'Value') returns toggle state of naivecheck


% --- Executes on button press in boostedcheck.
function boostedcheck_Callback(hObject, eventdata, handles)
% hObject    handle to boostedcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 clas11 = get(handles.boostedcheck,'Value');
if(clas11==1)
           
  end
% Hint: get(hObject,'Value') returns toggle state of boostedcheck


% --- Executes on button press in stumpscheck.
function stumpscheck_Callback(hObject, eventdata, handles)
% hObject    handle to stumpscheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clas10 = get(handles.stumpscheck,'Value');
if(clas10==1)
          
  end
% Hint: get(hObject,'Value') returns toggle state of stumpscheck


% --- Executes on button press in fccnncheck.
function fccnncheck_Callback(hObject, eventdata, handles)
% hObject    handle to fccnncheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 clas12 = get(handles.fccnncheck,'Value');
if(clas12==1)
         
  end
% Hint: get(hObject,'Value') returns toggle state of fccnncheck


% --- Executes on button press in qdacheck.
function qdacheck_Callback(hObject, eventdata, handles)
% hObject    handle to qdacheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 clas7 = get(handles.qdacheck,'Value');
if(clas7==1)
         
 end
% Hint: get(hObject,'Value') returns toggle state of qdacheck


% --- Executes on button press in rbfcheck.
function rbfcheck_Callback(hObject, eventdata, handles)
% hObject    handle to rbfcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clas8 = get(handles.rbfcheck,'Value');
if(clas8==1)
             
             hyperbf = inputdlg({'Lower Limit','Increment By','Upper Limit','Exponential Search'},...
              'RBF-SVM-C (Regularization paramter', [1 75; 1 75; 1 75; 1 75]);
handles.hyperbf=hyperbf;
guidata(hObject,handles);
assignin('base','hyperbf',handles.hyperbf);
hyperbf2 = inputdlg({'Lower Limit','Increment By','Upper Limit','Exponential Search'},...
              'RBF-SVM- Gamma', [1 75; 1 75; 1 75; 1 75]);
handles.hyperbf2=hyperbf2;
guidata(hObject,handles);
assignin('base','hyperbf2',handles.hyperbf2);
 end
% Hint: get(hObject,'Value') returns toggle state of rbfcheck


% --- Executes on button press in mlpcheck.
function mlpcheck_Callback(hObject, eventdata, handles)
% hObject    handle to mlpcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clas15 = get(handles.mlpcheck,'Value');
if(clas15==1)
            
   end
% Hint: get(hObject,'Value') returns toggle state of mlpcheck


% --- Executes on button press in rlrcheck.
function rlrcheck_Callback(hObject, eventdata, handles)
% hObject    handle to rlrcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clas17 = get(handles.rlrcheck,'Value');
if(clas17==1)
           
  end
% Hint: get(hObject,'Value') returns toggle state of rlrcheck


% --- Executes on button press in rvmcheck.
function rvmcheck_Callback(hObject, eventdata, handles)
% hObject    handle to rvmcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 clas19 = get(handles.rvmcheck,'Value');
if(clas19==1)
             
       hypervm = inputdlg({'Lower Limit','Increment By','Upper Limit','Exponential Search'},...
              'RVM- Gaussian kernel width', [1 75; 1 75; 1 75; 1 75]);
handles.hypervm=hypervm;
guidata(hObject,handles);
assignin('base','hypervm',handles.hypervm);
   end
% Hint: get(hObject,'Value') returns toggle state of rvmcheck


% --- Executes on button press in slrcheck.
function slrcheck_Callback(hObject, eventdata, handles)
% hObject    handle to slrcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 clas20 = get(handles.slrcheck,'Value');
if(clas20==1)
             
 end
% Hint: get(hObject,'Value') returns toggle state of slrcheck


% --- Executes on button press in consencheck.
function consencheck_Callback(hObject, eventdata, handles)
% hObject    handle to consencheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clas21 = get(handles.consencheck,'Value');
% Hint: get(hObject,'Value') returns toggle state of consencheck



function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
splirat = str2double(get( hObject, 'String' ));
handles.splirat=splirat;
guidata(hObject,handles);
assignin('base','splirat',handles.splirat);
% Hints: get(hObject,'String') returns contents of edit22 as text
%        str2double(get(hObject,'String')) returns contents of edit22 as a double


% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox34.
function checkbox34_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox.


% --- Executes on button press in checkbox55.
function checkbox55_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox55


% --- Executes on button press in checkbox57.
function checkbox57_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox57 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox57


% --- Executes on button press in checkbox58.
function checkbox58_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox58 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox58


% --- Executes on button press in checkbox59.
function checkbox59_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox59 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox59


% --- Executes on button press in checkbox60.
function checkbox60_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox60 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox60


% --- Executes on button press in checkbox61.
function checkbox61_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox61 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox61


% --- Executes on button press in checkbox62.
function checkbox62_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox62 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox62


% --- Executes on button press in checkbox63.
function checkbox63_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox63 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox63


% --- Executes on button press in checkbox64.
function checkbox64_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox64 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox64


% --- Executes on button press in checkbox65.
function checkbox65_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox65 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox65


% --- Executes on button press in checkbox66.
function checkbox66_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox66 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox66


% --- Executes on button press in checkbox67.
function checkbox67_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox67 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox67


% --- Executes on button press in checkbox68.
function checkbox68_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox68 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox68


% --- Executes on button press in checkbox69.
function checkbox69_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox69 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox69


% --- Executes on button press in checkbox70.
function checkbox70_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox70 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox70


% --- Executes on button press in checkbox71.
function checkbox71_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox71 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox71


% --- Executes on button press in checkbox72.
function checkbox72_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox72 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox72


% --- Executes on button press in checkbox73.
function checkbox73_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox73 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox73


% --- Executes on button press in checkbox74.
function checkbox74_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox74 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox74


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
folder_name=uigetdir('c:\');
%folder_name=strcat(savpthnme,savnme);
%[~,~,YourFile]=xlsread(test_Path);
handles.folder_name=folder_name;
guidata(hObject,handles);
assignin('base','folder_name',handles.folder_name);


% --- Executes on button press in checkbox75.
function checkbox75_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox75 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
exppara = get(handles.checkbox75,'Value');
% Hint: get(hObject,'Value') returns toggle state of checkbox75



function edit42_Callback(hObject, eventdata, handles)
% hObject    handle to edit42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lowvknn = str2double(get( hObject, 'String' ));
handles.lowvknn=lowvknn;
guidata(hObject,handles);
assignin('base','lowvknn',handles.lowvknn);
% Hints: get(hObject,'String') returns contents of edit42 as text
%        str2double(get(hObject,'String')) returns contents of edit42 as a double


% --- Executes during object creation, after setting all properties.
function edit42_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit43_Callback(hObject, eventdata, handles)
% hObject    handle to edit43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
incredknn = str2double(get( hObject, 'String' ));
handles.incredknn=incredknn;
guidata(hObject,handles);
assignin('base','incredknn',handles.incredknn);
% Hints: get(hObject,'String') returns contents of edit43 as text
%        str2double(get(hObject,'String')) returns contents of edit43 as a double


% --- Executes during object creation, after setting all properties.
function edit43_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit44_Callback(hObject, eventdata, handles)
% hObject    handle to edit44 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uppeknn = str2double(get( hObject, 'String' ));
handles.uppeknn=uppeknn;
guidata(hObject,handles);
assignin('base','uppeknn',handles.uppeknn);
% Hints: get(hObject,'String') returns contents of edit44 as text
%        str2double(get(hObject,'String')) returns contents of edit44 as a double


% --- Executes during object creation, after setting all properties.
function edit44_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit44 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox81.
function checkbox81_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox81 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox81


% --- Executes during object creation, after setting all properties.
function uibuttongroup36_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uibuttongroup36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in checkbox82.
function checkbox82_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox82 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox82


% --- Executes during object creation, after setting all properties.
function uibuttongroup7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uibuttongroup7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
