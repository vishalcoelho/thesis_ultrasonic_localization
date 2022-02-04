function varargout = gui_LocalizationServer(varargin)
%GUI_LOCALIZATIONSERVER M-file for gui_LocalizationServer.fig
%      GUI_LOCALIZATIONSERVER, by itself, creates a new GUI_LOCALIZATIONSERVER or raises the existing
%      singleton*.
%
%      H = GUI_LOCALIZATIONSERVER returns the handle to a new GUI_LOCALIZATIONSERVER or the handle to
%      the existing singleton*.
%
%      GUI_LOCALIZATIONSERVER('Property','Value',...) creates a new GUI_LOCALIZATIONSERVER using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to gui_LocalizationServer_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      GUI_LOCALIZATIONSERVER('CALLBACK') and GUI_LOCALIZATIONSERVER('CALLBACK',hObject,...) call the
%      local function named CALLBACK in GUI_LOCALIZATIONSERVER.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_LocalizationServer

% Last Modified by GUIDE v2.5 14-Nov-2009 12:27:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @gui_LocalizationServer_OpeningFcn, ...
    'gui_OutputFcn',  @gui_LocalizationServer_OutputFcn, ...
    'gui_LayoutFcn',  [], ...
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

% --- Executes just before gui_LocalizationServer is made visible.
function gui_LocalizationServer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

%Global directives
handles.DEBUG_MODE      = 0;
handles.SIMULATE        = 0;
handles.USE_POTENTIAL   = 1;
handles.USE_KALMAN      = 2;
handles.USE_METHOD      = handles.USE_POTENTIAL;

% Choose default command line output for gui_LocalizationServer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_LocalizationServer wait for user response (see UIRESUME)
% uiwait(handles.figure1);
CurrentSerialObjects = instrfindall; %Find all serial objects
if(isempty(CurrentSerialObjects)==0)
    fclose(CurrentSerialObjects); %Close all serial objects
    delete(CurrentSerialObjects); %Delete objects from memory
end

%Packet definition
handles.MaxNodes = 8;
handles.NumberOfNodes = get(handles.N_nodesMenu,'Value') + 2;
handles.NodesToTrack = get(handles.NTrackMenu,'Value');
handles.params_size = 17; %17 bytes of data per node
handles.DataPacketSize = (handles.params_size * handles.MaxNodes) + 4;
handles.HeaderSize = 4;     handles.TrailerSize = 2;
handles.TotalPacketSize = handles.HeaderSize + handles.DataPacketSize + handles.TrailerSize;

%Create the default serial port object
ser = serial(['COM' int2str(get(handles.ComPortMenu, 'Value'))], 'BaudRate', 115200, ...
    'OutputBufferSize', handles.TotalPacketSize * 20, 'InputBufferSize', handles.TotalPacketSize * 1000,...
    'Timeout', 1, 'FlowControl', 'none');
handles.AltSerialPort = ser;

%Configure a timer
% handles.PlotTimer = timer('BusyMode','drop','Period',20,...
% 'TimerFcn',{@PlotFcn,handles});

%Flags
StateMachine.Flags.ComPortOpen = 0;
StateMachine.Flags.StartLocalization = 0;
StateMachine.Flags.StartTracking = 0;
StateMachine.Flags.LocalizationDone = 0;
StateMachine.Flags.ComPortError = 0;

%Error Metrics
StateMachine.Stats.Points = 10;
StateMachine.Stats.ErrorMax = zeros([StateMachine.Stats.Points 1]); %Stores max relative error for nodes 3 to 8
StateMachine.Stats.ErrorAvg = zeros([StateMachine.Stats.Points 1]); %Stores avg relative error for nodes 3 to 8
StateMachine.Stats.Count = 1;
StateMachine.Stats.FilenameMax = ['MaxRelativeError',num2str(handles.NumberOfNodes),'.txt'];
StateMachine.Stats.FilenameAvg = ['AvgRelativeError',num2str(handles.NumberOfNodes),'.txt'];
StateMachine.Stats.fid1 = 0;
StateMachine.Stats.fid2 = 0;

%State Machine
StateMachine.State.Current = 1;
StateMachine.State.Next = 1;
StateMachine.Packet.Window = 10;
setappdata(handles.figure1,'StateMachine',StateMachine);

%Parameters
Matrices.Distance.Radio          = zeros(handles.MaxNodes,handles.MaxNodes,'uint16');
Matrices.Distance.Sound          = zeros(handles.MaxNodes,handles.MaxNodes,'uint16');
Matrices.TimeStamp.Radio         = zeros(handles.MaxNodes,handles.MaxNodes,'uint32');
Matrices.TimeStamp.Sound         = zeros(handles.MaxNodes,handles.MaxNodes,'uint32');
Matrices.RefNode.RSSI.Radio      = zeros(handles.MaxNodes,handles.MaxNodes,'uint8');
Matrices.RefNode.RSSI.Sync       = zeros(handles.MaxNodes,handles.MaxNodes,'uint8');
Matrices.RefNode.RSSI.BaseStation= zeros(handles.MaxNodes,handles.MaxNodes,'uint8');
Matrices.Mean.Radio              = zeros(handles.MaxNodes,handles.MaxNodes,'double');
Matrices.Mean.Sound              = zeros(handles.MaxNodes,handles.MaxNodes,'double');
Matrices.BaudRate                = [9600;19200;57600;115200;250000];
setappdata(hObject,'LocalizationData',Matrices);



%Save the handles
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = gui_LocalizationServer_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in StartLocalizationButton.
function StartLocalizationButton_Callback(hObject, eventdata, handles)
% hObject    handle to StartLocalizationButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Clear figure

clc;    %Clear the command line
StateMachine = getappdata(handles.figure1,'StateMachine');
if(StateMachine.Flags.ComPortOpen == 0)
    try
        fopen(handles.AltSerialPort);
    catch
        set(handles.StatusText,'String','Unable to Open COM Port');
        return;
    end
    StateMachine.Flags.ComPortOpen = 1;
end

%Reset Position vectors
Position.X = zeros([handles.MaxNodes , 1]);
Position.Y = zeros([handles.MaxNodes , 1]);
setappdata(handles.figure1,'NodeLocation',Position);

%Reset matrices
Matrices                         = getappdata(handles.figure1,'LocalizationData');
Matrices.Distance.Radio          = zeros(handles.MaxNodes,handles.MaxNodes,'uint16');
Matrices.Distance.Sound          = zeros(handles.MaxNodes,handles.MaxNodes,'uint16');
Matrices.TimeStamp.Radio         = zeros(handles.MaxNodes,handles.MaxNodes,'uint32');
Matrices.TimeStamp.Sound         = zeros(handles.MaxNodes,handles.MaxNodes,'uint32');
Matrices.RefNode.RSSI.Radio      = zeros(handles.MaxNodes,handles.MaxNodes,'uint8');
Matrices.RefNode.RSSI.Sync       = zeros(handles.MaxNodes,handles.MaxNodes,'uint8');
Matrices.RefNode.RSSI.BaseStation= zeros(handles.MaxNodes,handles.MaxNodes,'uint8');
Matrices.Mean.Radio              = zeros(handles.MaxNodes,handles.MaxNodes,'double');
Matrices.Mean.Sound              = zeros(handles.MaxNodes,handles.MaxNodes,'double');
setappdata(hObject,'LocalizationData',Matrices);

%Initialize debug data and packet count
if handles.DEBUG_MODE == 1
    StateMachine.DebugData = zeros([128 1]);
end
StateMachine.PacketCount = 0;
  
%Reset all associated flags, counters and states and save
StateMachine.Flags.StartLocalization = 0;
StateMachine.Flags.StartTracking = 0;
StateMachine.Flags.LocalizationDone = 0;
StateMachine.State.Current = 1;
StateMachine.State.Next =1;
StateMachine.Packet.Count = 1;
StateMachine.Stats.Count = 1;
setappdata(handles.figure1,'StateMachine',StateMachine);

try
    %RESET command to the master node
    reset = zeros(500,1,'uint8');
    fwrite(handles.AltSerialPort,reset,'uint8'); %Re-sync all nodes in the
    %event a node stopped middway during a prior transmission

    %STOP command to the master node
    header = zeros(4,1,'uint8');
    data = zeros(4,1,'uint8');
    trailer = zeros(2,1,'uint8');
    header(1) = hex2dec('7E');  %start byte
    header(2) = hex2dec('03');  %API computer command
    header(3) = 4;              %length LSB
    header(4) = 0;              %length MSB
    data(1) = hex2dec('01');    % command OPCODE - COMPUTER_COMMAND_OPCODE_SER_CTRL
    data(2) = hex2dec('00');    % transfer_ctrl
    trailer(1) = 255 - bitand(sum(data),255);   %crc
    trailer(2) = hex2dec('FF'); %end byte
    fwrite(handles.AltSerialPort,[header; data; trailer],'uint8');

    %FLUSH data in the receive buffer
    pause(0.2); %pause 200ms
    bytes_available = get(handles.AltSerialPort,'BytesAvailable');
    if bytes_available
        A = fread(handles.AltSerialPort,bytes_available,'uint8');
        %Reading the buffer clears it
    end

    %send START code to master node
    header(2) = hex2dec('03');  %API command
    header(3) = 4;              % length LSB
    header(4) = 0;              % length MSB
    data = zeros(4,1,'uint8');
    data(1) = hex2dec('01');    % command OPCODE - COMPUTER_COMMAND_OPCODE_SER_CTRL
    data(2) = hex2dec('01');    % transfer ctrl -> start transfers
    data(3) = uint8(handles.DataPacketSize);
    data(4) = 0;
    trailer(1) = 255 - bitand(sum(data),255); %crc check
    fwrite(handles.AltSerialPort,[header; data; trailer],'uint8');

    %Create an event each time a certain number of bytes becomes available
    %at the input buffer
    StateMachine.Flags.StartLocalization = 1;
catch
    StateMachine.Flags.StartLocalization = 0;
    set(handles.StatusText,'String','Error accessing the COM port.');
    setappdata(handles.figure1,'StateMachine',StateMachine);
    %fclose(handles.AltSerialPort);
    return;
end

%Close the COM port and set the BytesAvailable Fcn
if(StateMachine.Flags.ComPortOpen) %ComPort open?
    try
        fclose(handles.AltSerialPort);
        StateMachine.Flags.ComPortOpen = 0;
    catch
        set(handles.StatusText,'String','Localization failed: COM error');
    end
end

set(handles.AltSerialPort,'BytesAvailableFcnCount',handles.TotalPacketSize);
set(handles.AltSerialPort,'BytesAvailableFcnMode','byte')
set(handles.AltSerialPort,'BytesAvailableFcn',{@PopulateMatrix,handles});
guidata(hObject,handles);
try
    fopen(handles.AltSerialPort);
    StateMachine.Flags.ComPortOpen = 1;
catch
    set(handles.StatusText,'String','Unable to start localization');
    return;
end
setappdata(handles.figure1,'StateMachine',StateMachine);
guidata(hObject,handles);
set(handles.StatusText, 'String', 'Data written successfully, Localization initialized');


% --- Executes on button press in StartTrackingButton.
function StartTrackingButton_Callback(hObject, eventdata, handles)
% hObject    handle to StartTrackingButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% TRACKING THE TARGET
%Bounds check for number of trackable nodes

if (handles.NumberOfNodes - handles.NodesToTrack) <= 4
    set(handles.StatusText,'String','Not enough ground nodes');
    return;
end
switch(handles.USE_METHOD)
    case handles.USE_POTENTIAL
        StateMachine = getappdata(handles.figure1,'StateMachine');
        StateMachine.Flags.StartTracking = 1;
        setappdata(handles.figure1,'StateMachine',StateMachine);
    case handles.USE_KALMAN
        
end
guidata(hObject,handles);
%---END OF TRACKING---%

% --- Executes on selection change in N_nodesMenu.
function N_nodesMenu_Callback(hObject, eventdata, handles)
% hObject    handle to N_nodesMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = get(hObject,'String') returns N_nodesMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from N_nodesMenu
StateMachine = getappdata(handles.figure1,'StateMachine');
handles.NumberOfNodes = get(hObject,'Value')+2;
StateMachine.Stats.FilenameMax = ['MaxRelativeError',num2str(handles.NumberOfNodes),'.txt'];
StateMachine.Stats.FilenameAvg = ['AvgRelativeError',num2str(handles.NumberOfNodes),'.txt'];
try
    StateMachine.Stats.fid1 = fopen(StateMachine.Stats.FilenameMax,'w'); %open file for write
    StateMachine.Stats.fid2 = fopen(StateMachine.Stats.FilenameAvg,'w'); %open file for write
    setappdata(handles.figure1,'StateMachine',StateMachine);
    set(handles.StatusText,'String','Opened Error log successfully');
catch
    set(handles.StatusText,'String','Cant open error log');
    return;
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function N_nodesMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to N_nodesMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

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

% --- Executes on selection change in ComPortMenu.
function ComPortMenu_Callback(hObject, eventdata, handles)
% hObject    handle to ComPortMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = get(hObject,'String') returns ComPortMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ComPortMenu
handles.AvailableComPorts = getAvailableComPort;
set(hObject,'String',handles.AvailableComPorts);
%save the handles structure
guidata(hObject,handles);

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

% --- Executes during object creation, after setting all properties.
function DisplayAxes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DisplayAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: place code in OpeningFcn to populate DisplayAxes

function PopulateMatrix(hObject,eventdata,handles)
% hObject    Not Passed
% eventdata  Not Passed
% handles    handles to all UI objects and UI data

StateMachine = getappdata(handles.figure1,'StateMachine');

%Check for closed comport
if(StateMachine.Flags.ComPortOpen == 0)
    set(handles.StatusText,'String',['Com port closed, Counted ',...
        num2str(StateMachine.PacketCount),' packets']);
    return;
end

%Buffer stores 20 info packets and discards the rest, read 1 info packet at
%a time to make sure that u dont receive a packet mid transmission
BytesAvailable = get(handles.AltSerialPort,'BytesAvailable');
if(BytesAvailable)
    if(BytesAvailable >= handles.TotalPacketSize)
        %Read the first information packet
        A = uint8(fread(handles.AltSerialPort,handles.TotalPacketSize,'uint8'));
    else
        return;
    end
else
    return; %no data then return
end

%Check if localization was initialized
if(StateMachine.Flags.StartLocalization == 0) %Bogus data
    set(handles.StatusText,'String',['Data rejected, Localization complete ,Counted ',...
        num2str(StateMachine.PacketCount),' packets']);
    if(StateMachine.Flags.StartTracking == 0) %Dont track?
        return;
    end
end

try
    if (A(1) ~= hex2dec('7E')) || (A(2) ~= hex2dec('01')) || ((A(3) + A(4)*256) ~= handles.DataPacketSize) || (A(handles.TotalPacketSize) ~= hex2dec('FF'))
        set(handles.StatusText,'String','Error in the receieved packet form');
        fclose(handles.AltSerialPort);
        StateMachine.Flags.ComPortOpen = 0;
        %         StateMachine.Flags.ComPortError = 1;
        setappdata(handles.figure1,'StateMachine',StateMachine);
        %       guidata(hObject,handles);
        return;
    end
    if bitand(sum(A(5:(length(A)-2))), 255) + A(length(A)-1) ~= 255
        set(handles.StatusText, 'String', 'CRC error.');
        fclose(handles.AltSerialPort);
        StateMachine.Flags.ComPortOpen = 0;
        %         StateMachine.Flags.ComPortError = 1;
        setappdata(handles.figure1,'StateMachine',StateMachine);
        %       guidata(hObject,handles);
        return;
    end
catch
    set(handles.StatusText,'String','Error accessing the COM port');
    fclose(handles.AltSerialPort);
    StateMachine.Flags.ComPortOpen = 0;
    %     StateMachine.Flags.ComPortError = 1;
    setappdata(handles.figure1,'StateMachine',StateMachine);
    %   guidata(hObject,handles);
    return;
end
StateMachine.PacketCount = StateMachine.PacketCount + 1;
Matrices = getappdata(handles.figure1,'LocalizationData');
for i = 0 : handles.MaxNodes - 1 %hardware is zero indexed, 0 is 1 .... 7 is 8
    CurrentNode = i + 1;
    CurrentNodeOffset = i*handles.params_size;
    CurrentNodeStart = handles.HeaderSize + 1 + 4 + CurrentNodeOffset; %4 for the timestamp
    CurrentNodeEnd = CurrentNodeStart + handles.params_size - 1;
    CurrentNodeData = A(CurrentNodeStart : CurrentNodeEnd);
    params = CreateParams(CurrentNodeData);
    RefNodeRadio = params.reference_node_radio + 1; %Zero based in hardware
    if RefNodeRadio == 0
        RefNodeRadio = 1;
    end
    RefNodeSound = params.reference_node_sound + 1; %Zero based in hardware
    if RefNodeSound == 0
        RefNodeSound = 1;
    end
    %   Distances are in cms and time in us
    Matrices.Distance.Radio(CurrentNode,RefNodeRadio)            = params.radio_distance;
    Matrices.Distance.Radio(RefNodeRadio,CurrentNode)            = params.radio_distance;
    Matrices.Distance.Sound(CurrentNode,RefNodeSound)            = params.sound_distance;
    Matrices.Distance.Sound(RefNodeSound,CurrentNode)            = params.sound_distance;
    Matrices.TimeStamp.Radio(CurrentNode,RefNodeRadio)           = params.timestamp.radio;
    Matrices.TimeStamp.Sound(CurrentNode,RefNodeSound)           = params.timestamp.sound;
    Matrices.TimeStamp.Radio(RefNodeRadio,CurrentNode)           = params.timestamp.radio;
    Matrices.TimeStamp.Sound(RefNodeSound,CurrentNode)           = params.timestamp.sound;
    Matrices.RefNode.RSSI.Radio(CurrentNode,RefNodeRadio)        = params.rssi_reference_node_radio;
    Matrices.RefNode.RSSI.Sync(CurrentNode,RefNodeRadio)         = params.rssi_sync;
    Matrices.RefNode.RSSI.BaseStation(CurrentNode,RefNodeRadio)  = params.rssi_basestation;
    Matrices.RefNode.RSSI.Radio(RefNodeRadio,CurrentNode)        = params.rssi_reference_node_radio;
    Matrices.RefNode.RSSI.Sync(RefNodeRadio,CurrentNode)         = params.rssi_sync;
    Matrices.RefNode.RSSI.BaseStation(RefNodeRadio,CurrentNode)  = params.rssi_basestation;
end
if handles.SIMULATE == 1
    Matrices.RefNode.RSSI.BaseStation(1,1) = 65;
end
bar(handles.RSSIRadioAxes,Matrices.RefNode.RSSI.Radio(:,1));
bar(handles.RSSISyncAxes,Matrices.RefNode.RSSI.Sync(:,1));
bar(handles.RSSIBaseAxes,Matrices.RefNode.RSSI.BaseStation(:,1));
setappdata(handles.figure1,'LocalizationData',Matrices);
set(handles.StatusText, 'String', 'Data read successfully.');

%% Start of the State machine
% Matrices = getappdata(handles.figure1,'LocalizationData');
% StateMachine = getappdata(handles.figure1,'StateMachine');
%Average 10 packets to get stable data (assuming data is non zero)
if StateMachine.Flags.StartTracking == 0
    if (StateMachine.Packet.Count == 1)
        Matrices.Mean.Radio              = zeros(handles.MaxNodes,handles.MaxNodes,'double');
        Matrices.Mean.Sound              = zeros(handles.MaxNodes,handles.MaxNodes,'double');
    end
    if (StateMachine.Packet.Count <= StateMachine.Packet.Window)
        Matrices.Mean.Radio  = Matrices.Mean.Radio + double(Matrices.Distance.Radio);
        Matrices.Mean.Sound  = Matrices.Mean.Sound + double(Matrices.Distance.Sound);
        StateMachine.Packet.Count = StateMachine.Packet.Count + 1;
        setappdata(handles.figure1,'LocalizationData',Matrices);
        setappdata(handles.figure1,'StateMachine',StateMachine);
        return;
    else %get average data and reset the count
        Matrices.Mean.Radio  = Matrices.Mean.Radio./StateMachine.Packet.Window;
        Matrices.Mean.Sound  = Matrices.Mean.Sound./StateMachine.Packet.Window;
        StateMachine.Packet.Count = 1;
    end
else
    Matrices.Mean.Radio  = double(Matrices.Distance.Radio);
    Matrices.Mean.Sound  = double(Matrices.Distance.Sound);    
end
if(get(handles.RadioRadioButton,'Value') == get(handles.RadioRadioButton,'Max'))
    Rd = double(Matrices.Mean.Radio)./(2*100); %Work in meters and take into consideration RTT
else
    Rd = double(Matrices.Mean.Sound)./100;
end
tspan = [0 40];
if(StateMachine.Flags.LocalizationDone == 0) % Localization done ?
    switch(StateMachine.State.Current)
        case{1,2} %Dont run the dynamics, not enough data
            if Rd(1,2) ~= 0
                StateMachine.State.Next = StateMachine.State.Current + 1; %Next State
                %Displaying position data
                axes(handles.DisplayAxes)   %Assume the two nodes lie along a line
                plot(handles.DisplayAxes,[1,1],[Rd(1,2),1],'o');
            else
                StateMachine.State.Next = StateMachine.State.Current; %Stay in the same state
            end

        case 3
            if (length(nonzeros(Rd(3,1:3))) >= 2) % ?Sufficient data
                x0 = [rand(1); 0 ;rand(1); 0;...
                    rand(1); 0; rand(1); 0;...
                    rand(1); 0; rand(1); 0]; %Initial state
                %Prep the integrator for the dynamics
                [t,x] = ode15s(@(t,x) LocalizeNewNode(t,x,StateMachine.State.Current,Rd,handles),tspan,x0);
                [x_m x_n] = size(x);
                x_final = x(x_m,:); %Choose last row as final state
                for i = 1 : length(x_final)/4 %Save node positions
                    Position.X(i,1) = x_final(2*i-1 + 2*(i-1));
                    Position.Y(i,1) = x_final(2*i+1 + 2*(i-1));
                end
                setappdata(handles.figure1,'NodeLocation',Position);
                StateMachine.State.Next = StateMachine.State.Current + 1; %Next State
                %Displaying position data
                axes(handles.DisplayAxes)
                plot(handles.DisplayAxes,Position.X,Position.Y,'o')
            else    % ?Not enough data
                if handles.DEBUG_MODE == 1
                    %Keep a time record of two node measurements
                    temp = StateMachine.DebugData;
                    StateMachine.DebugData(1) = Rd(1,2);
                    StateMachine.DebugData(2:128) = temp(1:127);
                    %                     axes(handles.DebugAxes)
                    plot(handles.DebugAxes,[1:128],StateMachine.DebugData);

                    %Still display the last two nodes for debugging purposes
                    %Displaying position data
                    axes(handles.DisplayAxes)   %Assume the two nodes lie along a line
                    plot(handles.DisplayAxes,[1,1],[Rd(1,2),1],'o');
                end
                StateMachine.State.Next = StateMachine.State.Current; %stay in same state
            end

        case 4
            n = StateMachine.State.Current;
            if (n <= handles.NumberOfNodes) % ?Sufficient data
                if length(nonzeros(Rd(4,1:4))) >= 3 % ?Sufficient data
                    x0 = [rand(1); 0 ;rand(1); 0]; %Initial state
                    %Prep the integrator for the dynamics
                    [t,x] = ode15s(@(t,x) LocalizeNewNode(t,x,StateMachine.State.Current,Rd,handles),tspan,x0);
                    [x_m x_n] = size(x);
                    x_final = x(x_m,:); %Choose last row as final state
                    for i = 1 : length(x_final)/4 %Save node positions
                        X_node(i,1) = x_final(2*i-1 + 2*(i-1));
                        Y_node(i,1) = x_final(2*i+1 + 2*(i-1));
                    end
                    %Append position data and save
                    Position = getappdata(handles.figure1,'NodeLocation');
                    Position.X = [Position.X; X_node]; Position.Y = [Position.Y; Y_node];
                    setappdata(handles.figure1,'NodeLocation',Position);
                    StateMachine.State.Next = StateMachine.State.Current + 1; %Next State
                    %Displaying position data
                    axes(handles.DisplayAxes)
                    plot(handles.DisplayAxes,Position.X,Position.Y,'o');
                else
                    StateMachine.State.Next = StateMachine.State.Current; %stay in same state
                end
            else
                StateMachine.State.Next = 9;  %Go to final stage
            end

        case {5,6,7,8}
            n = StateMachine.State.Current;
            if (n <= handles.NumberOfNodes) % ?Sufficient data
                if length(nonzeros(Rd(n,1:n))) >= 3
                    x0 = [rand(1); 0 ;rand(1); 0]; %Initial state
                    %Prep the integrator for the dynamics
                    [t,x] = ode15s(@(t,x) LocalizeNewNode(t,x,StateMachine.State.Current,Rd,handles),tspan,x0);
                    [x_m x_n] = size(x);
                    x_final = x(x_m,:); %Choose last row as final state
                    for i = 1 : length(x_final)/4 %Save node positions
                        X_node(i,1) = x_final(2*i-1 + 2*(i-1));
                        Y_node(i,1) = x_final(2*i+1 + 2*(i-1));
                    end
                    %Append position data and save
                    Position = getappdata(handles.figure1,'NodeLocation');
                    Position.X = [Position.X; X_node]; Position.Y = [Position.Y; Y_node];
                    setappdata(handles.figure1,'NodeLocation',Position);
                    StateMachine.State.Next = StateMachine.State.Current + 1; %Next State
                    %Displaying position data
                    axes(handles.DisplayAxes)
                    plot(handles.DisplayAxes,Position.X,Position.Y,'o');

                else    % ?Not enough data
                    StateMachine.State.Next = StateMachine.State.Current; %stay in same state
                end
            else
                StateMachine.State.Next = 9;  %Go to final stage
            end

        case 9 %all nodes
            Position = getappdata(handles.figure1,'NodeLocation');
            x0 = zeros([handles.NumberOfNodes*4 1]);
            for i = 1 : length(Position.X) %Initial States
                x0(2*i-1 + 2*(i-1)) = Position.X(i);  x0(2*i+1 + 2*(i-1)) = Position.Y(i);
                x0(2*i + 2*(i-1))   = 0;              x0(2*i + 2 + 2*(i-1))   = 0;
            end
            %Prep the integrator for the dynamics
            [t,x] = ode15s(@(t,x) LocalizeNewNode(t,x,StateMachine.State.Current,Rd,handles),tspan,x0);
            [x_m x_n] = size(x);
            x_final = x(x_m,:); %Choose last row as final state
            for i = 1 : length(x_final)/4 %Save node positions
                X_node(i,1) = x_final(2*i-1 + 2*(i-1));
                Y_node(i,1) = x_final(2*i+1 + 2*(i-1));
            end
            %Rewrite the finalised data
            Position.X = X_node; Position.Y = Y_node;
            setappdata(handles.figure1,'NodeLocation',Position);

            %Displaying position data
            axes(handles.DisplayAxes)
            plot(handles.DisplayAxes,Position.X,Position.Y,'o');

            %States and flags
            StateMachine.State.Next = 1; %Revert back to initial state
            StateMachine.Flags.LocalizationDone = 1; %Localization complete
            StateMachine.Flags.StartLocalization = 0; %push start button to reactivate localization

            %Completion Message
            set(handles.StatusText,'String','Localization Completed');
            
            if (handles.DEBUG_MODE == 1)
                if(StateMachine.State.Current > 2)
                    %Error metrics
                    Position = getappdata(handles.figure1,'NodeLocation');
                    for i = 1 : length(Position.X)
                        for k = 1 : length(Position.X)
                            Rd_hat(i,k) = sqrt((Position.X(i)-Position.X(k))^2 + (Position.Y(i)-Position.Y(k))^2); %cm's
                            Error(i,k) = abs(Rd(i,k) - Rd_hat(i,k));
                            if(Rd(i,k) ~= 0)
                                RelativeError(i,k) = (Error(i,k)/Rd(i,k))*100;
                            else
                                RelativeError(i,k) = 0;
                            end
                        end
                    end
                    clc
                    disp('Rd'),disp(Rd)
                    disp('Rd_hat'),disp(Rd_hat)
                    disp('Error'),disp(Error)
                    disp('Relative Error'),disp(RelativeError)
                    %store max relative error
                    row = StateMachine.Stats.Count;
                    StateMachine.Stats.ErrorMax(row,1) = max(max(RelativeError));
                    StateMachine.Stats.ErrorAvg(row,1) = mean(mean(RelativeError));
                    StateMachine.Stats.Count = StateMachine.Stats.Count+1;
                    if StateMachine.Stats.Count >= StateMachine.Stats.Points
                        try
                            StateMachine.Stats.Count = 1; %reset
                            setappdata(handles.figure1,'StateMachine',StateMachine);
                            dlmwrite(StateMachine.Stats.FilenameMax,StateMachine.Stats.ErrorMax,'');
                            dlmwrite(StateMachine.Stats.FilenameAvg,StateMachine.Stats.ErrorAvg,'');
                            fclose(StateMachine.Stats.FilenameMax);
                            fclose(StateMachine.Stats.FilenameAvg);
                            set(handles.StatusText,'String','Wrote error log successfully');
                            disp('!!!Wrote error log successfully!!!');
                         catch
                            set(handles.StatusText,'String','Cant write error log');
                        end
                    end
                end
            end
        otherwise %Do Nothing
    end
    StateMachine.State.Current = StateMachine.State.Next;
    setappdata(handles.figure1,'StateMachine',StateMachine);

    %--------------EOSM------------%
    % End of the state machine
end
if(StateMachine.Flags.StartTracking == 1)
    handles.NodesToTrack = get(handles.NTrackMenu,'Value');
    Position = getappdata(handles.figure1,'NodeLocation');
    n = handles.NumberOfNodes -handles.NodesToTrack + 1;
    for i = 1 : handles.NodesToTrack %Initial States
        x0 = zeros([4 1]);
        x0 = [Position.X((i-1) + n ); 0; Position.Y((i-1) + n); 0];
       
        %Prep the integrator for the dynamics
        [t,x] = ode15s(@(t,x) LocalizeNewNode(t,x,n+i-1,Rd,handles),tspan,x0);
        [x_m x_n] = size(x);
        x_final = x(x_m,:); %Choose last row as final state
        for j = 1 : length(x_final)/4 %Save node positions
            X_node(j,1) = x_final(2*j-1 + 2*(j-1));
            Y_node(j,1) = x_final(2*j+1 + 2*(j-1));
        end
        %Rewrite the finalised data
        Position.X(n+i-1) = X_node; Position.Y(n+i-1) = Y_node;
        setappdata(handles.figure1,'NodeLocation',Position);
    end
    %Displaying position data
    axes(handles.DisplayAxes)
    plot(handles.DisplayAxes,Position.X,Position.Y,'ro');
end
% if strcmp(get(handles.PlotTimer,'Running'),'off')
%     start(handles.PlotTimer);
% end
function params = CreateParams( A )
A = uint8(A);
params.radio_distance              = typecast(A(1:2), 'uint16');
params.sound_distance              = typecast(A(3:4), 'uint16');
params.timestamp.radio             = typecast(A(5:8), 'uint32');
params.timestamp.sound             = typecast(A(9:12), 'uint32');
params.reference_node_radio        = typecast(A(13), 'uint8');
params.reference_node_sound        = typecast(A(14), 'uint8');
params.rssi_reference_node_radio   = typecast(A(15), 'uint8');
params.rssi_sync                   = typecast(A(16), 'uint8');
params.rssi_basestation            = typecast(A(17), 'uint8');


function xdot = LocalizeNewNode(t,x,State,Rd,handles)
%t time span for integration
%x State Vector
%State state of the machine
%Rd radio or ultrasound distances

switch(State)
    case {1,2} %Wait three states for atleast 3 nodes information to arrive
        xdot = x; %Keep initial states, this case should never be executed
    case 3 % 1st 3 nodes can be localized
        %Renaming position vectors
        x1 = x(1);  y1 = x(3);
        x2 = x(5);  y2 = x(7);
        x3 = x(9);  y3 = x(11);
        %Renaming velocity vectors
        x1dot=x(2);     y1dot=x(4);
        x2dot=x(6);     y2dot=x(8);
        x3dot=x(10);    y3dot=x(12);
        % Calculation of distance  between Nodes
        % R21=R12; if there is no measurement noise(N), In this case N=0
        r12=sqrt((x1-x2)^2+(y1-y2)^2);      r21=r12;
        r13=sqrt((x1-x3)^2+(y1-y3)^2);      r31=r13;
        r23=sqrt((x2-x3)^2+(y2-y3)^2);      r32=r23;
        %Gains
        k = 1;          kv = 1;
        % Calculating force for Node 1 in X & Y direction
        F1x= k*(r12-Rd(1,2))*((x2-x1)/r12) - kv*x1dot + k*(r13-Rd(1,3))*((x3-x1)/r13);
        F1y= k*(r12-Rd(1,2))*((y2-y1)/r12) - kv*y1dot + k*(r13-Rd(1,3))*((y3-y1)/r13);

        % Calculating force for Node 2 X & Y direction
        F2x= k*(r21-Rd(2,1))*((x1-x2)/r21) - kv*x2dot + k*(r23-Rd(2,3))*((x3-x2)/r23);
        F2y= k*(r21-Rd(2,1))*((y1-y2)/r21) - kv*y2dot + k*(r23-Rd(2,3))*((y3-y2)/r23);

        % Calculating force for Node 3 X & Y direction
        F3x= k*(r31-Rd(3,1))*((x1-x3)/r31) - kv*x3dot + k*(r32-Rd(3,2))*((x2-x3)/r32);
        F3y= k*(r31-Rd(3,1))*((y1-y3)/r31) - kv*y3dot + k*(r32-Rd(3,2))*((y2-y3)/r32);

        % Sensor Node Dynamics (Xdot=Ax+Bu)
        xdot=[x(2);  % x1
            F1x;
            x(4);  % y1
            F1y;   % End of First Node Dynamics
            x(6);  % x2
            F2x;
            x(8);  % y2
            F2y;   % End of Second Node Dynamics
            x(10); % x3
            F3x;
            x(12); % y3
            F3y];   %% End of Third Node Dynamics

    case 4 %4th node
        %The final state vector from the previous
        %state gives the final position of the nodes
        %under consideration in that state(state machine)
        %Retrieve position of previously localized nodes
        Position = getappdata(handles.figure1,'NodeLocation');
        X_f = Position.X;   Y_f = Position.Y;
        %Renaming the positions
        x4 = x(1);  y4 = x(3);
        %Renaming the velocities
        x4dot = x(2); y4dot = x(4);
        %Calculation of distance between Nodes
        r14=sqrt((X_f(1)-x4)^2+(Y_f(1)-y4)^2);              r41=r14;
        r12=sqrt((X_f(1)-X_f(2))^2+(Y_f(1)-Y_f(2))^2);      r21=r12;
        r13=sqrt((X_f(1)-X_f(3))^2+(Y_f(1)-Y_f(3))^2);      r31=r13;
        r24=sqrt((X_f(2)-x4)^2+(Y_f(2)-y4)^2);              r42=r24;
        r34=sqrt((X_f(3)-x4)^2+(Y_f(3)-y4)^2);              r43=r34;

        a= r14^2+r24^2-r12^2;
        a1= r14^2+r34^2-r13^2;

        b= 2*r14*r24;
        b1= 2*r14*r34;
        theta42d=acos(a/b);
        theta43d=acos(a1/b1);
        cosTHETA42=a/b;
        cosTHETA43=a1/b1;
        delta_abx= (2*b*((x4-X_f(1))+(x4-X_f(2))) - (2*a/(r14*r24))*(r24^2*(x4-X_f(1)) + r14^2*(x4-X_f(2))))/(4*r14^2*r24^2);
        delta_aby= (2*b*((y4-Y_f(1)+(y4-Y_f(2))) - (2*a/(r14*r24))*(r24^2*(y4-Y_f(1)) + r14^2*(y4-Y_f(2)))))/(4*r14^2*r24^2);

        delta_abx1= (2*b1*((x4-X_f(1))+(x4-X_f(3))) - (2*a1/(r14*r34))*(r34^2*(x4-X_f(1)) + r14^2*(x4-X_f(3))))/(4*r14^2*r34^2);
        delta_aby1= (2*b1*((y4-Y_f(1))+(y4-Y_f(3))) - (2*a1/(r14*r34))*(r34^2*(y4-Y_f(1)) + r14^2*(y4-Y_f(3))))/(4*r14^2*r34^2);
        c=a/b;
        c1=a1/b1;
        d=-1/(sqrt(1-c^2));
        d1=-1/(sqrt(1-c1^2));
        % Gains
        k=1;
        % Dissipative gain
        kv=1;

        % Calculating force for Node 1 in X & Y direction
        F4x = k*(r41-Rd(4,1))*((X_f(1)-x4)/r41) + k*(r42-Rd(4,2))*((X_f(2)-x4)/r42) + k*(r43-Rd(4,3))*((X_f(3)-x4)/r43) - 1*1*(cosTHETA42-cos(theta42d))*delta_abx - 1*1*(cosTHETA43-cos(theta43d))*delta_abx1 - kv*x4dot;
        F4y = k*(r41-Rd(4,1))*((Y_f(1)-y4)/r41) + k*(r42-Rd(4,2))*((Y_f(2)-y4)/r42) + k*(r43-Rd(4,3))*((Y_f(3)-y4)/r43) - 1*1*(cosTHETA42-cos(theta42d))*delta_aby - 1*1*(cosTHETA43-cos(theta43d))*delta_aby1 - kv*y4dot;

        % Sensor Node Dynamics (Xdot=Ax+Bu)
        xdot=[ x(2);  % x1
            F4x;
            x(4);  % y1
            F4y];   % End of First Node Dynamics

    case {5,6,7,8} %5th 6th 7th and 8th node and end state
        n = State;
        %The final state vector from the previous
        %state gives the final position of the nodes
        %under consideration in that state(state machine)
        %Retrieve position of previously localized nodes
        Position = getappdata(handles.figure1,'NodeLocation');
        X_f = Position.X;   Y_f = Position.Y;
        %Renaming position variables for the nth node
        xn = x(1);      yn = x(3);
        %Renaming velocity variables
        xndot = x(2);   yndot = x(4);
        Fnx = 0; Fny = 0;
        % Gains
        k=1;            kv=1;
        %         count = 1;
        %         if count <= 3
        for i = 1:(n-1)
            if(Rd(n,i) ~= 0)
                R(n,i) = sqrt((xn-X_f(i))^2 + (yn-Y_f(i))^2);
                Fnx = Fnx + k*(R(n,i)-Rd(n,i))*((X_f(i)-xn)/R(n,i));
                Fny = Fny + k*(R(n,i)-Rd(n,i))*((Y_f(i)-yn)/R(n,i));
                %                     count = count + 1;
            end
        end
        %         end
        Fnx = Fnx - kv*xndot;
        Fny = Fny - kv*yndot;
        xdot=[ x(2);  % x1
            Fnx;
            x(4);  % y1
            Fny];

    case 9
        % All nodes
        n = handles.NumberOfNodes;
        %The final state vector from the previous
        %state gives the final position of the nodes
        %under consideration in that state(state machine)
        StateVector = getappdata(handles.figure1,'StateVector');
        %Retrieve position of previously localized nodes
        Position = getappdata(handles.figure1,'NodeLocation');
        X_f = Position.X;   Y_f = Position.Y;
        xdot = zeros([handles.NumberOfNodes*4 1]);
        R    = zeros([handles.NumberOfNodes,handles.NumberOfNodes]);
        for q = 1 : n
            %Renaming position variables for the nth node
            xq = x(2*q-1 + 2*(q-1));      yq = x(2*q+1 + 2*(q-1));
            %Renaming velocity variables
            xqdot = x(2*q + 2*(q-1));   yqdot = x(2*q + 2 + 2*(q-1));
            Fqx = 0; Fqy = 0;
            % Gains
            k=1;            kv=1;
            for i = 1:n %Each node will have n-1 neighbours i ~= q
                %                 count = 1;
                if i ~= q
                    %                     if count <= 3 %Calculate forces for only 3 neighbours
                    if (Rd(q,i) ~= 0) % Only nonzero data
                        R(q,i) = sqrt((xq-X_f(i))^2 + (yq-Y_f(i))^2);
                        Fqx = Fqx + k*(R(q,i)-Rd(q,i))*((X_f(i)-xq)/R(q,i));
                        Fqy = Fqy + k*(R(q,i)-Rd(q,i))*((Y_f(i)-yq)/R(q,i));
                        %                             count = count + 1;
                    end
                    %                     end
                end
            end
            Fqx = Fqx - kv*xqdot;
            Fqy = Fqy - kv*yqdot;
            xdot(2*q-1 + 2*(q-1)) = xqdot; xdot(2*q+1 + 2*(q-1))= yqdot;
            xdot(2*q + 2*(q-1)) = Fqx; xdot(2*q + 2 + 2*(q-1)) = Fqy;
        end

    otherwise %Default state;

end

% --- Executes on button press in DisconnectButton.
function ConnectButton_Callback(hObject, eventdata, handles)
% hObject    handle to ConnectButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(handles.AltSerialPort == 0)
    hndAltSerialPort = serial(['COM' int2str(get(handles.ComPortMenu,'Value'))],'BaudRate', 115200, ...
        'OutputBufferSize', 2048, 'InputBufferSize', 2048, 'Timeout', 1, 'FlowControl', 'none');
    set(hndAltSerialPort ,'BytesAvailableFcnCount',1);
    set(hndAltSerialPort ,'BytesAvailableFcnMode','byte')
    set(hndAltSerialPort ,'BytesAvailableFcn',{@PopulateMatrix,handles});
    handles.AltSerialPort = hndAltSerialPort;
end
%Flush the buffer
BytesAvailable = get(handles.AltSerialPort,'BytesAvailable');
if(BytesAvailable)
    A = fread(handles.AltSerialPort,BytesAvailable);
end
Matrices = getappdata(handles.figure1,'LocalizationData');
StateMachine = getappdata(handles.figure1,'StateMachine');
BaudRate = Matrices.BaudRate;
ChosenComPort = handles.AvailableComPorts(get(handles.ComPortMenu,'Value'));
set(handles.AltSerialPort,'Port',cell2mat(ChosenComPort));
set(handles.AltSerialPort,'BaudRate',BaudRate(get(handles.BaudRateMenu,'Value')));
guidata(hObject,handles);
try
    fopen(handles.AltSerialPort);
    StateMachine.Flags.ComPortOpen = 1;
    set(handles.StatusText,'String',['Successfully opened ',...
        get(handles.AltSerialPort,'Port'),' at ',...
        int2str(get(handles.AltSerialPort,'BaudRate')),' baud']);
catch
    set(handles.StatusText,'String',['Unable to open ',...
        get(handles.AltSerialPort,'Port'),' at ',...
        int2str(get(handles.AltSerialPort,'BaudRate')),' baud']);
    return;
end
%Save the state machine
setappdata(handles.figure1,'StateMachine',StateMachine);
%Save the handles structure
guidata(hObject,handles);

% --- Executes on button press in DisconnectButton.
function DisconnectButton_Callback(hObject, eventdata, handles)
% hObject    handle to DisconnectButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
StateMachine = getappdata(handles.figure1,'StateMachine');
if(StateMachine.Flags.ComPortOpen == 0)
    set(handles.StatusText,'String','Port already disconnected');
    return;
end
try
    bytes_available = get(handles.AltSerialPort,'BytesAvailable');
    if bytes_available
        A = fread(handles.AltSerialPort,bytes_available,'uint8');
        %Reading the buffer clears it
    end
    fclose(handles.AltSerialPort); %Close the port
    delete(handles.AltSerialPort); %delete serial object
    handles.AltSerialPort = 0;
    StateMachine.Flags.ComPortOpen = 0; %set com to 0
    handles.AvailableComPorts = getAvailableComPort; %Re populate com menu
    set(handles.ComPortMenu,'String',handles.AvailableComPorts);
    set(handles.StatusText,'String','Port closed successfully');
catch
    set(handles.StatusText,'String','Unable to close port');
    return;
end

%save handles structure
guidata(hObject,handles);

% --- Executes on button press in RadioRadioButton.
function RadioRadioButton_Callback(hObject, eventdata, handles)
% hObject    handle to RadioRadioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RadioRadioButton
if(get(hObject,'Value') == get(hObject,'Max'))
    set(handles.SoundRadioButton,'Value',get(handles.SoundRadioButton,'Min'));
end

% --- Executes on button press in SoundRadioButton.
function SoundRadioButton_Callback(hObject, eventdata, handles)
% hObject    handle to SoundRadioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SoundRadioButton
if(get(hObject,'Value') == get(hObject,'Max'))
    set(handles.RadioRadioButton,'Value',get(handles.RadioRadioButton,'Min'));
end

% --- Executes on button press in ShowNodesButton.
function ShowNodesButton_Callback(hObject, eventdata, handles)
% hObject    handle to ShowNodesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
StateMachine = getappdata(handles.figure1,'StateMachine');
if StateMachine.Flags.LocalizationDone == 0
    set(handles.StatusText,'String','Localization not complete, displaying available data');
    axes(handles.DisplayAxes)
    %Get Positional data
    Position = getappdata(handles.figure1,'NodeLocation');
    for i = 1 : handles.MaxNodes %No need to plot again
        if (Position.X(i) ~= 0 && Position.Y(i) ~= 0)
            grid on;
            text(Position.X(i),Position.Y(i),num2str(i-1),....
                'HorizontalAlignment','left','VerticalAlignment','top');
        end
    end
else %Localization is complete
    axes(handles.DisplayAxes)
    %Get Positional data
    Position = getappdata(handles.figure1,'NodeLocation');
    for i = 1 : handles.NumberOfNodes %No need to plot again
        text(Position.X(i),Position.Y(i),num2str(i-1),....
            'HorizontalAlignment','left','VerticalAlignment','top');
    end
end


% --- Executes on selection change in NTrackMenu.
function NTrackMenu_Callback(hObject, eventdata, handles)
% hObject    handle to NTrackMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns NTrackMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from NTrackMenu


% --- Executes during object creation, after setting all properties.
function NTrackMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NTrackMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%---Function to reset the serial port---%
function status = SerialPortReset(hndSerialPort,handles)
%hndSerialPort - handle to the serial port object
%handles - GUI data
%status is 1 if reset is successful and 0 otherwise

status = 0;
try
    fclose(hndSerialPort);
    set(handles.StatusText,'String','Port closed successfully');
catch
    set(handles.StatusText,'String','Unable to close port');
    status = 0;
    return;
end
%Once closed set the bytes available function
set(hndSerialPort,'BytesAvailableFcnCount',handles.TotalPacketSize);
set(hndSerialPort,'BytesAvailableFcnMode','byte')
set(hndSerialPort,'BytesAvailableFcn',{@PopulateMatrix,handles});
%Reopen the port
try
    fopen(hndSerialPort);
    set(handles.StatusText,'String','Reopened COM port');
    %RESET command to the master node
    reset = zeros(500,1,'uint8');
    fwrite(hndSerialPort,reset,'uint8'); %Re-sync all nodes in the
    %event a node stopped middway during a prior transmission

    %STOP command to the master node
    header = zeros(4,1,'uint8');
    data = zeros(4,1,'uint8');
    trailer = zeros(2,1,'uint8');
    header(1) = hex2dec('7E');  %start byte
    header(2) = hex2dec('03');  %API computer command
    header(3) = 4;              %length LSB
    header(4) = 0;              %length MSB
    data(1) = hex2dec('01');    % command OPCODE - COMPUTER_COMMAND_OPCODE_SER_CTRL
    data(2) = hex2dec('00');    % transfer_ctrl
    trailer(1) = 255 - bitand(sum(data),255);   %crc
    trailer(2) = hex2dec('FF'); %end byte
    fwrite(hndSerialPort,[header; data; trailer],'uint8');

    %FLUSH data in the receive buffer
    pause(0.2); %pause 200ms
    bytes_available = get(hndSerialPort,'BytesAvailable');
    if bytes_available
        A = fread(hndSerialPort,bytes_available,'uint8');
        %Reading the buffer clears it
    end

    %send START code to master node
    header(2) = hex2dec('03');  %API command
    header(3) = 4;              % length LSB
    header(4) = 0;              % length MSB
    data = zeros(4,1,'uint8');
    data(1) = hex2dec('01');    % command OPCODE - COMPUTER_COMMAND_OPCODE_SER_CTRL
    data(2) = hex2dec('01');    % transfer ctrl -> start transfers
    data(3) = uint8(handles.DataPacketSize);
    data(4) = 0;
    trailer(1) = 255 - bitand(sum(data),255); %crc check
    fwrite(hndSerialPort,[header; data; trailer],'uint8');
    status = 1;
catch
    set(handles.StatusText,'String','Failed to reopen COM port');
    status = 0;
    return;
end

function PlotFcn(hObject, eventdata, handles)
% hObject    handle to the display axes(see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
Position = getappdata(handles.figure1,'NodeLocation');
plot(Position.X,Position.Y,'o');
hold on;
for i = 1 : handles.MaxNodes %No need to plot again
        if (Position.X(i) ~= 0 && Position.Y(i) ~= 0)
            grid on;
            text(Position.X(i),Position.Y(i),num2str(i),....
                'HorizontalAlignment','left','VerticalAlignment','top');
        end
    end
hold off;