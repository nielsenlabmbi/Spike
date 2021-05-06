function varargout = sortGui(varargin)
% SORTGUI MATLAB code for sortGui.fig
%      SORTGUI, by itself, creates a new SORTGUI or raises the existing
%      singleton*.
%
%      H = SORTGUI returns the handle to a new SORTGUI or the handle to
%      the existing singleton*.
%
%      SORTGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SORTGUI.M with the given input arguments.
%
%      SORTGUI('Property','Value',...) creates a new SORTGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before sortGui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sortGui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sortGui

% Last Modified by GUIDE v2.5 25-Jan-2019 16:07:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sortGui_OpeningFcn, ...
                   'gui_OutputFcn',  @sortGui_OutputFcn, ...
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


% --- Executes just before sortGui is made visible.
function sortGui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sortGui (see VARARGIN)

% Choose default command line output for sortGui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes sortGui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = sortGui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes during object creation, after setting all properties.
function unitsMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to unitsMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function propOneMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to propOneMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function propTwoMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to propTwoMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function compUnitMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to compUnitMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in loadFile.
function loadFile_Callback(hObject, eventdata, handles)
global Properties idk Time CHsPos PathName FileName
load([fileparts(mfilename('fullpath')) '/Settings.mat'])
[FileName,PathName,FilterIndex] = uigetfile;
load([PathName FileName])

CHsPos = 14;
Amp0 = 5;
Time = 15;
for i = 1:length(PropTitles)
  switch PropTitles{i}
      case 'CHs pos'
          CHsPos=i;
      case 'Amp (0)'
          Amp0=i;
      case 'Time'
          Time = i;
  end
end
idk = sortTheUnits(Properties,idk,CHsPos);

handles.propOneMenu.String = PropTitles;
handles.propOneMenu.String = PropTitles;
handles.fileTxt.String = FileName;
handles.propTwoMenu.String = PropTitles;
handles.propTwoMenu.Value = 5;
units = sort(unique(idk));
handles.unitsMenu.String = units;
handles.unitsMenu.Value = 1;
handles.compUnitMenu.String = units;
handles.compUnitMenu.Value = 1;

meanFR = (sum(idk==units(handles.unitsMenu.Value))./max(Properties(Time,:)))*30000;
handles.meanFRAns.String = num2str(meanFR);

Tst = sort(Properties(Time,idk==units(1)));
dTst = diff(Tst);
cont = sum(dTst<36)./sum(idk==units(1));
handles.ISIContAns.String = num2str(cont);

scatter(handles.theOnlyAxes,Properties(CHsPos,(idk==units(1))),Properties(Amp0,(idk==units(1))),'SizeData',1,'CData',[0 0 0],'Marker','.')
% hObject    handle to loadFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in prevUnit.
function prevUnit_Callback(hObject, eventdata, handles)
global Properties idk Time
units = sort(unique(idk));
currVal = handles.unitsMenu.Value;
currVal = currVal-1;
if currVal<1
    currVal = length(handles.unitsMenu.String);
end
handles.unitsMenu.Value = currVal;

meanFR = (sum(idk==units(handles.unitsMenu.Value))./max(Properties(Time,:)))*30000;
handles.meanFRAns.String = num2str(meanFR);

Tst = sort(Properties(Time,idk==units(handles.unitsMenu.Value)));
dTst = diff(Tst);
cont = sum(dTst<36)./sum(idk==units(handles.unitsMenu.Value));
handles.ISIContAns.String = num2str(round(cont,3,'significant'));
scatter(handles.theOnlyAxes,Properties(handles.propOneMenu.Value,(idk==units(handles.unitsMenu.Value))),Properties(handles.propTwoMenu.Value,(idk==units(handles.unitsMenu.Value))),'SizeData',1,'CData',[0 0 0],'Marker','.')
% hObject    handle to prevUnit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in nextUnit.
function nextUnit_Callback(hObject, eventdata, handles)
global Properties idk Time
units = sort(unique(idk));
currVal = handles.unitsMenu.Value;
currVal = currVal+1;
if currVal>length(handles.unitsMenu.String)
    currVal = 1;
end
handles.unitsMenu.Value = currVal;

meanFR = (sum(idk==units(handles.unitsMenu.Value))./max(Properties(Time,:)))*30000;
handles.meanFRAns.String = num2str(meanFR);

Tst = sort(Properties(Time,idk==units(handles.unitsMenu.Value)));
dTst = diff(Tst);
cont = sum(dTst<36)./sum(idk==units(handles.unitsMenu.Value));
handles.ISIContAns.String = num2str(round(cont,3,'significant'));

scatter(handles.theOnlyAxes,Properties(handles.propOneMenu.Value,(idk==units(handles.unitsMenu.Value))),Properties(handles.propTwoMenu.Value,(idk==units(handles.unitsMenu.Value))),'SizeData',1,'CData',[0 0 0],'Marker','.')
% hObject    handle to nextUnit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on selection change in unitsMenu.
function unitsMenu_Callback(hObject, eventdata, handles)
global Properties idk Time
units = sort(unique(idk));

meanFR = (sum(idk==units(handles.unitsMenu.Value))./max(Properties(Time,:)))*30000;
handles.meanFRAns.String = num2str(meanFR);

Tst = sort(Properties(Time,idk==units(handles.unitsMenu.Value)));
dTst = diff(Tst);
cont = sum(dTst<36)./sum(idk==units(handles.unitsMenu.Value));
handles.ISIContAns.String = num2str(round(cont,3,'significant'));

scatter(handles.theOnlyAxes,Properties(handles.propOneMenu.Value,(idk==units(handles.unitsMenu.Value))),Properties(handles.propTwoMenu.Value,(idk==units(handles.unitsMenu.Value))),'SizeData',1,'CData',[0 0 0],'Marker','.')
% hObject    handle to unitsMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns unitsMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from unitsMenu


% --- Executes on button press in nextPropOne.
function nextPropOne_Callback(hObject, eventdata, handles)
global Properties idk
units = sort(unique(idk));
currVal = handles.propOneMenu.Value;
currVal = currVal+1;
if currVal>length(handles.propOneMenu.String)
    currVal = 1;
end
handles.propOneMenu.Value = currVal;
scatter(handles.theOnlyAxes,Properties(handles.propOneMenu.Value,(idk==units(handles.unitsMenu.Value))),Properties(handles.propTwoMenu.Value,(idk==units(handles.unitsMenu.Value))),'SizeData',1,'CData',[0 0 0],'Marker','.')
% hObject    handle to nextPropOne (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in prevPropOne.
function prevPropOne_Callback(hObject, eventdata, handles)
global Properties idk
units = sort(unique(idk));
currVal = handles.propOneMenu.Value;
currVal = currVal-1;
if currVal<1
    currVal = length(handles.propOneMenu.String);
end
handles.propOneMenu.Value = currVal;
scatter(handles.theOnlyAxes,Properties(handles.propOneMenu.Value,(idk==units(handles.unitsMenu.Value))),Properties(handles.propTwoMenu.Value,(idk==units(handles.unitsMenu.Value))),'SizeData',1,'CData',[0 0 0],'Marker','.')
% hObject    handle to prevPropOne (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in propOneMenu.
function propOneMenu_Callback(hObject, eventdata, handles)
global Properties idk
units = sort(unique(idk));
scatter(handles.theOnlyAxes,Properties(handles.propOneMenu.Value,(idk==units(handles.unitsMenu.Value))),Properties(handles.propTwoMenu.Value,(idk==units(handles.unitsMenu.Value))),'SizeData',1,'CData',[0 0 0],'Marker','.')
% hObject    handle to propOneMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns propOneMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from propOneMenu


% --- Executes on button press in nextPropTwo.
function nextPropTwo_Callback(hObject, eventdata, handles)
global Properties idk
units = sort(unique(idk));
currVal = handles.propTwoMenu.Value;
currVal = currVal+1;
if currVal>length(handles.propTwoMenu.String)
    currVal = 1;
end
handles.propTwoMenu.Value = currVal;
scatter(handles.theOnlyAxes,Properties(handles.propOneMenu.Value,(idk==units(handles.unitsMenu.Value))),Properties(handles.propTwoMenu.Value,(idk==units(handles.unitsMenu.Value))),'SizeData',1,'CData',[0 0 0],'Marker','.')

% hObject    handle to nextPropTwo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in prevPropTwo.
function prevPropTwo_Callback(hObject, eventdata, handles)
global Properties idk
units = sort(unique(idk));
currVal = handles.propTwoMenu.Value;
currVal = currVal-1;
if currVal<1
    currVal = length(handles.propTwoMenu.String);
end
handles.propTwoMenu.Value = currVal;
scatter(handles.theOnlyAxes,Properties(handles.propOneMenu.Value,(idk==units(handles.unitsMenu.Value))),Properties(handles.propTwoMenu.Value,(idk==units(handles.unitsMenu.Value))),'SizeData',1,'CData',[0 0 0],'Marker','.')
% hObject    handle to prevPropTwo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in propTwoMenu.
function propTwoMenu_Callback(hObject, eventdata, handles)
global Properties idk
units = sort(unique(idk));
scatter(handles.theOnlyAxes,Properties(handles.propOneMenu.Value,(idk==units(handles.unitsMenu.Value))),Properties(handles.propTwoMenu.Value,(idk==units(handles.unitsMenu.Value))),'SizeData',1,'CData',[0 0 0],'Marker','.')
% hObject    handle to propTwoMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns propTwoMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from propTwoMenu


% --- Executes on button press in uinitISI.
function uinitISI_Callback(hObject, eventdata, handles)
global Properties idk Time
units = sort(unique(idk));
Tst = sort(Properties(Time,idk==units(handles.unitsMenu.Value)));
dTst = diff(Tst);
hist(handles.theOnlyAxes,(dTst(dTst<120))./30,24);
hold on
plot(handles.theOnlyAxes,[1.2 1.2],[0 max(hist((dTst(dTst<120))./30,24))]);
hold off
set(handles.theOnlyAxes,'XLim',[0 4])
% hObject    handle to uinitISI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in cleanUnit.
function cleanUnit_Callback(hObject, eventdata, handles)
global Properties idk Time
units = sort(unique(idk));
scatter(handles.theOnlyAxes,Properties(handles.propOneMenu.Value,(idk==units(handles.unitsMenu.Value))),Properties(handles.propTwoMenu.Value,(idk==units(handles.unitsMenu.Value))),'SizeData',1,'CData',[0 0 0],'Marker','.')
BW = impoly;
nodes=getPosition(BW);
selected = (inpolygon(Properties(handles.propOneMenu.Value,idk==units(handles.unitsMenu.Value)),Properties(handles.propTwoMenu.Value,idk==units(handles.unitsMenu.Value)),nodes(:,1),nodes(:,2)));
Ploted = find(idk==units(handles.unitsMenu.Value));
idk(Ploted(~selected))=0;

if sum(idk==units(handles.unitsMenu.Value)) == 0
for unit2move = units(handles.unitsMenu.Value)+1:max(idk) 
    idk(idk==unit2move) = unit2move-1;
end
end

meanFR = (sum(idk==units(handles.unitsMenu.Value))./max(Properties(Time,:)))*30000;
handles.meanFRAns.String = num2str(meanFR);

Tst = sort(Properties(Time,idk==units(handles.unitsMenu.Value)));
dTst = diff(Tst);
cont = sum(dTst<36)./sum(idk==units(handles.unitsMenu.Value));
handles.ISIContAns.String = num2str(round(cont,3,'significant'));
scatter(handles.theOnlyAxes,Properties(handles.propOneMenu.Value,(idk==units(handles.unitsMenu.Value))),Properties(handles.propTwoMenu.Value,(idk==units(handles.unitsMenu.Value))),'SizeData',1,'CData',[0 0 0],'Marker','.')
% hObject    handle to cleanUnit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in sortUnit.
function sortUnit_Callback(hObject, eventdata, handles)
global Properties idk CHsPos Time
units = sort(unique(idk));
scatter(handles.theOnlyAxes,Properties(handles.propOneMenu.Value,(idk==units(handles.unitsMenu.Value))),Properties(handles.propTwoMenu.Value,(idk==units(handles.unitsMenu.Value))),'SizeData',1,'CData',[0 0 0],'Marker','.')
BW = impoly;
nodes=getPosition(BW);
selected = (inpolygon(Properties(handles.propOneMenu.Value,idk==units(handles.unitsMenu.Value)),Properties(handles.propTwoMenu.Value,idk==units(handles.unitsMenu.Value)),nodes(:,1),nodes(:,2)));
Ploted = find(idk==units(handles.unitsMenu.Value));
idk(Ploted(selected))=max(idk)+1;

if sum(idk==units(handles.unitsMenu.Value)) == 0
for unit2move = units(handles.unitsMenu.Value)+1:max(idk) 
    idk(idk==unit2move) = unit2move-1;
end
end

newunit = max(idk);
idkpos = find(idk==newunit,1);
idk = sortTheUnits(Properties,idk,CHsPos);
newunit = idk(idkpos);

units = sort(unique(idk));
handles.unitsMenu.String = units;
handles.compUnitMenu.String = units;
handles.unitsMenu.Value = newunit+1;

meanFR = (sum(idk==units(handles.unitsMenu.Value))./max(Properties(Time,:)))*30000;
handles.meanFRAns.String = num2str(meanFR);

Tst = sort(Properties(Time,idk==units(handles.unitsMenu.Value)));
dTst = diff(Tst);
cont = sum(dTst<36)./sum(idk==units(handles.unitsMenu.Value));
handles.ISIContAns.String = num2str(round(cont,3,'significant'));
scatter(handles.theOnlyAxes,Properties(handles.propOneMenu.Value,(idk==units(handles.unitsMenu.Value))),Properties(handles.propTwoMenu.Value,(idk==units(handles.unitsMenu.Value))),'SizeData',1,'CData',[0 0 0],'Marker','.')

% hObject    handle to sortUnit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in delteUnit.
function delteUnit_Callback(hObject, eventdata, handles)
global Properties idk CHsPos Time

units = sort(unique(idk));
idk(idk==units(handles.unitsMenu.Value))=0;

for unit2move = units(handles.unitsMenu.Value)+1:max(idk) 
    idk(idk==unit2move) = unit2move-1;
end

units = sort(unique(idk));
handles.unitsMenu.String = units;
handles.unitsMenu.Value = 1;
handles.compUnitMenu.String = units;
handles.compUnitMenu.Value = 1;

meanFR = (sum(idk==units(handles.unitsMenu.Value))./max(Properties(Time,:)))*30000;
handles.meanFRAns.String = num2str(meanFR);

Tst = sort(Properties(Time,idk==units(handles.unitsMenu.Value)));
dTst = diff(Tst);
cont = sum(dTst<36)./sum(idk==units(handles.unitsMenu.Value));
handles.ISIContAns.String = num2str(round(cont,3,'significant'));
scatter(handles.theOnlyAxes,Properties(handles.propOneMenu.Value,(idk==units(handles.unitsMenu.Value))),Properties(handles.propTwoMenu.Value,(idk==units(handles.unitsMenu.Value))),'SizeData',1,'CData',[0 0 0],'Marker','.')

% hObject    handle to delteUnit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in nextCompUnit.
function nextCompUnit_Callback(hObject, eventdata, handles)
global Properties idk
units = sort(unique(idk));
currVal = handles.compUnitMenu.Value;
currVal = currVal+1;
if currVal>length(handles.compUnitMenu.String)
    currVal = 1;
end
handles.compUnitMenu.Value = currVal;

units = sort(unique(idk));
if ~handles.allUnitRadio.Value  
hold on
scatter(handles.theOnlyAxes,Properties(handles.propOneMenu.Value,(idk==units(handles.compUnitMenu.Value))),Properties(handles.propTwoMenu.Value,(idk==units(handles.compUnitMenu.Value))),'SizeData',1,'CData',[1 0 0],'Marker','.')
hold off
end
% hObject    handle to nextCompUnit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in prevCompUnit.
function prevCompUnit_Callback(hObject, eventdata, handles)
global Properties idk
units = sort(unique(idk));
currVal = handles.compUnitMenu.Value;
currVal = currVal-1;
if currVal<1
    currVal = length(handles.compUnitMenu.String);
end
handles.compUnitMenu.Value = currVal;

units = sort(unique(idk));
if ~handles.allUnitRadio.Value  
hold on
scatter(handles.theOnlyAxes,Properties(handles.propOneMenu.Value,(idk==units(handles.compUnitMenu.Value))),Properties(handles.propTwoMenu.Value,(idk==units(handles.compUnitMenu.Value))),'SizeData',1,'CData',[1 0 0],'Marker','.')
hold off
end
% hObject    handle to prevCompUnit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on selection change in compUnitMenu.
function compUnitMenu_Callback(hObject, eventdata, handles)
global Properties idk
units = sort(unique(idk));
if ~handles.allUnitRadio.Value  
hold on
scatter(handles.theOnlyAxes,Properties(handles.propOneMenu.Value,(idk==units(handles.compUnitMenu.Value))),Properties(handles.propTwoMenu.Value,(idk==units(handles.compUnitMenu.Value))),'SizeData',1,'CData',[1 0 0],'Marker','.')
hold off
end
% hObject    handle to compUnitMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns compUnitMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from compUnitMenu


% --- Executes on button press in compUnitPlot.
function compUnitPlot_Callback(hObject, eventdata, handles)
global Properties idk CHsPos Time
units = sort(unique(idk));
if handles.allUnitRadio.Value  
hold on
for unit = 1:length(units)
if unit ~= handles.unitsMenu.Value
scatter(handles.theOnlyAxes,Properties(handles.propOneMenu.Value,(idk==units(unit))),Properties(handles.propTwoMenu.Value,(idk==units(unit))),'SizeData',1,'Marker','.')
end
end
hold off
else
hold on
scatter(handles.theOnlyAxes,Properties(handles.propOneMenu.Value,(idk==units(handles.compUnitMenu.Value))),Properties(handles.propTwoMenu.Value,(idk==units(handles.compUnitMenu.Value))),'SizeData',1,'CData',[1 0 0],'Marker','.')
hold off
end
% hObject    handle to compUnitPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in compUnitCrossISI.
function compUnitCrossISI_Callback(hObject, eventdata, handles)
global Properties idk CHsPos Time
units = sort(unique(idk));
unit1 = units(handles.unitsMenu.Value);
unit2 = units(handles.compUnitMenu.Value);
unit1Tst = sort(Properties(Time,idk==unit1));
unit2Tst = sort(Properties(Time,idk==unit2));
Tst = [unit1Tst unit2Tst];
Source = [zeros(1,length(unit1Tst)) ones(1,length(unit2Tst),1)];
[Tst Order] = sort(Tst);
Source = Source(Order);
dTst = diff(Tst);
dSource = diff(Source);
dTst(dSource == 0)=[];
hist(handles.theOnlyAxes,(dTst(dTst<240))./30,24);
hold on
plot(handles.theOnlyAxes,[1.2 1.2],[0 max(hist((dTst(dTst<120))./30,24))]);
hold off
% hObject    handle to compUnitCrossISI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 

% --- Executes on button press in compUnitMerge.
function compUnitMerge_Callback(hObject, eventdata, handles)
global Properties idk CHsPos Time

units = sort(unique(idk));
idk(idk==units(handles.unitsMenu.Value))=units(handles.compUnitMenu.Value);

for unit2move = units(handles.unitsMenu.Value)+1:max(idk) 
    idk(idk==unit2move) = unit2move-1;
end

units = sort(unique(idk));
handles.unitsMenu.String = units;
handles.unitsMenu.Value = 1;
handles.compUnitMenu.String = units;
if handles.compUnitMenu.Value<handles.unitsMenu.Value
handles.unitsMenu.Value = handles.compUnitMenu.Value;
else
handles.unitsMenu.Value = handles.compUnitMenu.Value;
end

meanFR = (sum(idk==units(handles.unitsMenu.Value))./max(Properties(Time,:)))*30000;
handles.meanFRAns.String = num2str(meanFR);

Tst = sort(Properties(Time,idk==units(handles.unitsMenu.Value)));
dTst = diff(Tst);
cont = sum(dTst<36)./sum(idk==units(handles.unitsMenu.Value));
handles.ISIContAns.String = num2str(round(cont,3,'significant'));
scatter(handles.theOnlyAxes,Properties(handles.propOneMenu.Value,(idk==units(handles.unitsMenu.Value))),Properties(handles.propTwoMenu.Value,(idk==units(handles.unitsMenu.Value))),'SizeData',1,'CData',[0 0 0],'Marker','.')
% hObject    handle to compUnitMerge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in saveFile.
function saveFile_Callback(hObject, eventdata, handles)
global Properties idk Time CHsPos PathName FileName
save([PathName FileName],'idk','Properties','-append')
% hObject    handle to saveFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function allUnitRadio_Callback(hObject, eventdata, handles, varargin)


% --- Executes during object creation, after setting all properties.
function theOnlyAxes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to theOnlyAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate theOnlyAxes
