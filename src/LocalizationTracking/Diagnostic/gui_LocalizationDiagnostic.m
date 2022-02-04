function varargout = gui_LocalizationDiagnostic(varargin)
% GUI_LOCALIZATIONDIAGNOSTIC M-file for gui_LocalizationDiagnostic.fig
%      GUI_LOCALIZATIONDIAGNOSTIC, by itself, creates a new GUI_LOCALIZATIONDIAGNOSTIC or raises the existing
%      singleton*.
%
%      H = GUI_LOCALIZATIONDIAGNOSTIC returns the handle to a new GUI_LOCALIZATIONDIAGNOSTIC or the handle to
%      the existing singleton*.
%
%      GUI_LOCALIZATIONDIAGNOSTIC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_LOCALIZATIONDIAGNOSTIC.M with the given input arguments.
%
%      GUI_LOCALIZATIONDIAGNOSTIC('Property','Value',...) creates a new GUI_LOCALIZATIONDIAGNOSTIC or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_LocalizationDiagnostic_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_LocalizationDiagnostic_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_LocalizationDiagnostic

% Last Modified by GUIDE v2.5 01-Nov-2009 15:17:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_LocalizationDiagnostic_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_LocalizationDiagnostic_OutputFcn, ...
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


% --- Executes just before gui_LocalizationDiagnostic is made visible.
function gui_LocalizationDiagnostic_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_LocalizationDiagnostic (see VARARGIN)

% Choose default command line output for gui_LocalizationDiagnostic
handles.output = hObject;
CurrentSerialObjects = instrfindall; %Find all serial objects
if(isempty(CurrentSerialObjects)==0)
    fclose(CurrentSerialObjects); %Close all serial objects
    delete(CurrentSerialObjects); %Delete objects from memory
end
hndSerialPort = serial(['COM' int2str(get(handles.ComPortMenu,'Value'))],'BaudRate', 115200, ...
    'OutputBufferSize', 2048, 'InputBufferSize', 2048, 'Timeout', 1, 'FlowControl', 'none');
set(hndSerialPort ,'BytesAvailableFcnCount',1);
set(hndSerialPort ,'BytesAvailableFcnMode','byte')
set(hndSerialPort ,'BytesAvailableFcn',{@SerialReceive,handles});
handles.SerialPort               = hndSerialPort;
handles.Flags.ComPortOpen        = 0;
handles.MaxNodes                 = 8;
handles.NumberOfNodes            = get(handles.NNodesMenu,'Value') + 2;
handles.ParametersSize           = 17; %17 bytes of data per node
handles.DataPacketSize           = uint16(handles.ParametersSize * handles.MaxNodes);
Matrices.Distance.Radio          = zeros(handles.MaxNodes,handles.MaxNodes,'uint16');
Matrices.Distance.Sound          = zeros(handles.MaxNodes,handles.MaxNodes,'uint16');
Matrices.TimeStamp.Radio         = zeros(handles.MaxNodes,handles.MaxNodes,'uint32');
Matrices.TimeStamp.Sound         = zeros(handles.MaxNodes,handles.MaxNodes,'uint32');
Matrices.RefNode.RSSI.Radio      = zeros(handles.MaxNodes,handles.MaxNodes,'uint8');
Matrices.RefNode.RSSI.Sync       = zeros(handles.MaxNodes,handles.MaxNodes,'uint8');
Matrices.RefNode.RSSI.BaseStation= zeros(handles.MaxNodes,handles.MaxNodes,'uint8');
Matrices.BaudRate                = [9600;19200;57600;115200;250000];
setappdata(hObject,'LocalizationData',Matrices);
% Delta step
handles.Step = 1;
%String index
handles.Count = 1;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_LocalizationDiagnostic wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_LocalizationDiagnostic_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in NNodesMenu.
function NNodesMenu_Callback(hObject, eventdata, handles)
% hObject    handle to NNodesMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns NNodesMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from NNodesMenu
handles.NumberOfNodes = get(hObject,'Value') + 2;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function NNodesMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NNodesMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in TransmitPacketButton.
function TransmitPacketButton_Callback(hObject, eventdata, handles)
% hObject    handle to TransmitPacketButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
A = CreatePacket(handles);
handles.Step = handles.Step + 1;
if(handles.Step > handles.NumberOfNodes - 1)
    handles.Step = 1;
end
%Save the handles structure
guidata(hObject,handles);
try
    fwrite(handles.SerialPort,A,'uint8');
    set(handles.StatusText,'String','Packet successfully transmitted');
catch
    if(handles.Flags.ComPortOpen == 0)
        set(handles.StatusText,'String','COM port not open');
    else
        set(handles.StatusText,'String','Unable to transmit packet');
    end
    return;
end
% --- Executes on selection change in ComPortMenu.
function ComPortMenu_Callback(hObject, eventdata, handles)
% hObject    handle to ComPortMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns ComPortMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ComPortMenu


% --- Executes during object creation, after setting all properties.
function ComPortMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ComPortMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.AvailableComPorts = getAvailableComPort;
set(hObject,'String',handles.AvailableComPorts);
%save the handles structure
guidata(hObject,handles);
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in BaudRateMenu.
function BaudRateMenu_Callback(hObject, eventdata, handles)
% hObject    handle to BaudRateMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns BaudRateMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from BaudRateMenu


% --- Executes during object creation, after setting all properties.
function BaudRateMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BaudRateMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ConnectButton.
function ConnectButton_Callback(hObject, eventdata, handles)
% hObject    handle to ConnectButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(handles.SerialPort == 0)
    hndSerialPort = serial(['COM' int2str(get(handles.ComPortMenu,'Value'))],'BaudRate', 115200, ...
    'OutputBufferSize', 2048, 'InputBufferSize', 2048, 'Timeout', 1, 'FlowControl', 'none');
    set(hndSerialPort ,'BytesAvailableFcnCount',1);
    set(hndSerialPort ,'BytesAvailableFcnMode','byte')
    set(hndSerialPort ,'BytesAvailableFcn',{@SerialReceive,handles});
    handles.SerialPort = hndSerialPort
end
try
    Matrices = getappdata(handles.figure1,'LocalizationData');
    BaudRate = Matrices.BaudRate;
    ChosenComPort = handles.AvailableComPorts(get(handles.ComPortMenu,'Value'));
    set(handles.SerialPort,'Port',cell2mat(ChosenComPort));
    set(handles.SerialPort,'BaudRate',BaudRate(get(handles.BaudRateMenu,'Value')));
    guidata(hObject,handles);
    fopen(handles.SerialPort);
    handles.Flags.ComPortOpen = 1;
    set(handles.StatusText,'String',['Successfully opened ',...
        get(handles.SerialPort,'Port'),' at ',...
        int2str(get(handles.SerialPort,'BaudRate')),' baud']);
catch
    set(handles.StatusText,'String',['Unable to open ',...
        get(handles.SerialPort,'Port'),' at ',...
        int2str(get(handles.SerialPort,'BaudRate')),' baud']);
    return;    
end

%Save the handles structure
guidata(hObject,handles);


% --- Executes on selection change in ViewDataMenu.
function ViewDataMenu_Callback(hObject, eventdata, handles)
% hObject    handle to ViewDataMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns ViewDataMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ViewDataMenu
Matrices = getappdata(handles.figure1,'LocalizationData');
index = get(hObject,'Value');
switch index
    case 1      
        set(handles.ViewDataTable,'Data',Matrices.Distance.Radio);
    case 2      
        set(handles.ViewDataTable,'Data',Matrices.Distance.Sound);
    case 3      
        set(handles.ViewDataTable,'Data',Matrices.TimeStamp.Radio);
    case 4      
        set(handles.ViewDataTable,'Data',Matrices.TimeStamp.Sound);
    case 5      
        set(handles.ViewDataTable,'Data',Matrices.RefNode.RSSI.Radio);
    case 6      
        set(handles.ViewDataTable,'Data',Matrices.RefNode.RSSI.Sync);
    case 7      
        set(handles.ViewDataTable,'Data',Matrices.RefNode.RSSI.BaseStation);
    otherwise
        set(handles.ViewDataTable,'Data',Matrices.Distance.Radio);
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function ViewDataMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ViewDataMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function ViewDataTable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ViewDataTable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when entered data in editable cell(s) in ViewDataTable.
function ViewDataTable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to ViewDataTable (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
index = get(handles.ViewDataMenu,'Value');
Matrices = getappdata(handles.figure1,'LocalizationData');
Data = get(hObject,'Data');
switch(index)
    case 1
        Matrices.Distance.Radio          = Data;
    case 2
        Matrices.Distance.Sound          = Data;
    case 3
        Matrices.TimeStamp.Radio         = Data;
    case 4
        Matrices.TimeStamp.Sound         = Data;
    case 5
        Matrices.RefNode.RSSI.Radio      = Data;
    case 6
        Matrices.RefNode.RSSI.Sync       = Data;
    case 7
        Matrices.RefNode.RSSI.BaseStation= Data;
    otherwise
end
setappdata(handles.figure1,'LocalizationData',Matrices);
guidata(hObject,handles);

function Packet = CreatePacket(handles)
% hObject    Not used
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure that has the handles to all GUI objects
Step = handles.Step;
Matrices = getappdata(handles.figure1,'LocalizationData');
Packet = zeros(handles.DataPacketSize + 6,1,'uint8');
Packet(1)  = hex2dec('7E');  %Start byte
Packet(2)  = hex2dec('01');  %api command
DataPacketSize = typecast(handles.DataPacketSize,'uint8');
Packet(3:4) = DataPacketSize;
for i = 1 : handles.NumberOfNodes
   CurNode = i;
   RefNode = CurNode - Step;
   if(RefNode <= 0) 
       RefNode = RefNode + handles.NumberOfNodes;
   end
   Start = handles.ParametersSize * (i - 1) + 5;
   RadioDistance   = uint16(Matrices.Distance.Radio(RefNode,CurNode));%Radio Distance
   SoundDistance   = uint16(Matrices.Distance.Sound(RefNode,CurNode));%Sound Distance
   TimeStampRadio  = uint32(Matrices.TimeStamp.Radio(RefNode,CurNode));%Radio timestamp
   TimeStampSound  = uint32(Matrices.TimeStamp.Sound(RefNode,CurNode));%Sound timestamp
   RefNodeSound    = uint8(RefNode);%Ref node sound
   RefNodeRadio    = uint8(RefNode);%Ref node radio
   RSSIRadio       = uint8(Matrices.RefNode.RSSI.Radio(RefNode,CurNode));%Ref node RSSI
   RSSISync        = uint8(Matrices.RefNode.RSSI.Sync(RefNode,CurNode));%Ref sync RSSI
   RSSIBase        = uint8(Matrices.RefNode.RSSI.BaseStation(RefNode,CurNode));%Ref basestation RSSI
   Packet(Start:Start+1)   = typecast(RadioDistance,'uint8');
   Packet(Start+2:Start+3) = typecast(SoundDistance ,'uint8');
   Packet(Start+4:Start+7) = typecast(TimeStampRadio ,'uint8');
   Packet(Start+8:Start+11)= typecast(TimeStampSound ,'uint8');
   Packet(Start+12)        = typecast(RefNodeSound ,'uint8');
   Packet(Start+13)        = typecast(RefNodeRadio ,'uint8');
   Packet(Start+14)        = typecast(RSSIRadio ,'uint8');
   Packet(Start+15)        = typecast(RSSISync ,'uint8');
   Packet(Start+16)        = typecast(RSSIBase ,'uint8');
end
crcsum = sum(Packet(5:handles.DataPacketSize+5));
Packet(handles.DataPacketSize + 5) = 255 - bitand(crcsum,255); %crc
Packet(handles.DataPacketSize + 6) = hex2dec('FF'); %Packet end



% --- Executes on button press in GenerateDataButton.
function GenerateDataButton_Callback(hObject, eventdata, handles)
% hObject    handle to GenerateDataButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Reset Step information
handles.Step = 1;
Matrices = getappdata(handles.figure1,'LocalizationData');
Matrices.Distance.Radio          = zeros(handles.MaxNodes,handles.MaxNodes,'uint16');
Matrices.Distance.Sound          = zeros(handles.MaxNodes,handles.MaxNodes,'uint16');
Matrices.TimeStamp.Radio         = zeros(handles.MaxNodes,handles.MaxNodes,'uint32');
Matrices.TimeStamp.Sound         = zeros(handles.MaxNodes,handles.MaxNodes,'uint32');
Matrices.RefNode.RSSI.Radio      = zeros(handles.MaxNodes,handles.MaxNodes,'uint8');
Matrices.RefNode.RSSI.Sync       = zeros(handles.MaxNodes,handles.MaxNodes,'uint8');
Matrices.RefNode.RSSI.BaseStation= zeros(handles.MaxNodes,handles.MaxNodes,'uint8');
Data = Matrices.Distance.Radio;
for i = 1 : handles.NumberOfNodes
    for j = i : handles.NumberOfNodes
        if(i == j)
            Data(i,j) = 0;
        else
            Data(i,j) = ceil(rand(1) * 1000);
            Data(j,i) = Data(i,j);
        end      
    end
end
Matrices.Distance.Radio = Data;
set(handles.ViewDataTable,'Data',Matrices.Distance.Radio);
%Save the handles structure
guidata(hObject,handles);
%Save application data
setappdata(handles.figure1,'LocalizationData',Matrices);


% --- Executes on button press in TestButton.
function TestButton_Callback(hObject, eventdata, handles)
% hObject    handle to TestButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(handles.Flags.ComPortOpen == 0)
    set(handles.StatusText,'String','COM port not open');
    return;
end
handles.Count = handles.Count + 1;
if handles.Count > 6
    handles.Count = 1;
end
TestString = {'Localization Diagnostic v1.1';
              'Hello World!';
              'Creature Feature';
              'Taxidermy';
              'Prolonged exposure';
              'Kalman Filter'};
try
    fwrite(handles.SerialPort,char(TestString(handles.Count)),'uchar');
    set(handles.StatusText,'String','Test string successfully sent');
catch
    set(handles.StatusText,'String','Unable to send test string');
    return;
end
    
    
%save the handles structure
guidata(hObject,handles);

function SerialReceive(hObject, eventdata, handles)
% hObject    none (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% %Disable the BytesAvailableFcn
% set(handles.SerialPort,'BytesAvailableFcnCount',0);
% set(handles.SerialPort,'BytesAvailableFcn',0);
% guidata(hObject,handles);
% 
% %Read the buffer
% A = fread(handles.SerialPort,'uint8');
% set(handles.StatusText,'String',int2str(A));
% %save the handles structure
% guidata(hObject,handles);


% --- Executes on button press in DisconnectButton.
function DisconnectButton_Callback(hObject, eventdata, handles)
% hObject    handle to DisconnectButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(handles.Flags.ComPortOpen == 0)
    set(handles.StatusText,'String','Port already disconnected');
    return;
end
try
    bytes_available = get(handles.SerialPort,'BytesAvailable');
    if bytes_available
        A = fread(handles.SerialPort,bytes_available,'uint8');
        %Reading the buffer clears it
    end
    fclose(handles.SerialPort); %Close the port
    delete(handles.SerialPort); %delete serial object
    handles.SerialPort = 0;
    handles.Flags.ComPortOpen = 0; %set com to 0
    handles.AvailableComPorts = getAvailableComPort; %Re populate com menu
    set(handles.ComPortMenu,'String',handles.AvailableComPorts);
    set(handles.StatusText,'String','Port closed successfully');
catch
   set(handles.StatusText,'String','Unable to close port');
   return;    
end

%save handles structure
guidata(hObject,handles);

