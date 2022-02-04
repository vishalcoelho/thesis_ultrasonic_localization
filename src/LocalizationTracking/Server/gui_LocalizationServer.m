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

% Last Modified by GUIDE v2.5 01-Nov-2009 16:21:19

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

%Create the default serial port object
ser = serial(['COM' int2str(get(handles.ComPortMenu, 'Value'))], 'BaudRate', 115200, ...
    'OutputBufferSize', 2048, 'InputBufferSize', 2048, 'Timeout', 1, 'FlowControl', 'none');
handles.AltSerialPort = ser;
handles.MaxNodes = 8;
handles.NumberOfNodes = get(handles.N_nodesMenu,'Value') + 2;
handles.params_size = 17; %17 bytes of data per node
handles.DataPacketSize = handles.params_size * handles.MaxNodes;
handles.Flags.ComPortOpen = 0;
% set(handles.AltSerialPort,'BytesAvailableFcnCount',handles.DataPacketSize+6);
% set(handles.AltSerialPort,'BytesAvailableFcnMode','byte')
% set(handles.AltSerialPort,'BytesAvailableFcn',{@PopulateMatrix,handles});
handles.Flags.StartLocalization = 0;
handles.Flags.LocalizationDone = 0;
StateMachine.State.Current = 1;
StateMachine.State.Next = 1;
setappdata(handles.figure1,'StateMachine',StateMachine);
Matrices.Distance.Radio          = zeros(handles.MaxNodes,handles.MaxNodes,'uint16');
Matrices.Distance.Sound          = zeros(handles.MaxNodes,handles.MaxNodes,'uint16');
Matrices.TimeStamp.Radio         = zeros(handles.MaxNodes,handles.MaxNodes,'uint32');
Matrices.TimeStamp.Sound         = zeros(handles.MaxNodes,handles.MaxNodes,'uint32');
Matrices.RefNode.RSSI.Radio      = zeros(handles.MaxNodes,handles.MaxNodes,'uint8');
Matrices.RefNode.RSSI.Sync       = zeros(handles.MaxNodes,handles.MaxNodes,'uint8');
Matrices.RefNode.RSSI.BaseStation= zeros(handles.MaxNodes,handles.MaxNodes,'uint8');
Matrices.BaudRate                = [9600;19200;57600;115200;250000];
setappdata(hObject,'LocalizationData',Matrices);
guidata(hObject, handles);
%Localization Data
% Initial Points - We need atleast 3 nodes to start a localization
% algotrithm and the assumption is we have this data at the start
Vector.State.Initial = zeros(handles.NumberOfNodes*2,2,'uint16');
Vector.State.Final   = zeros(handles.NumberOfNodes*2,2,'uint16');
x_i(1) =  ceil(rand(1)*1000); y_i(1) =  ceil(rand(1)*1000);
x_i(2) =  ceil(rand(1)*1000); y_i(2) =  ceil(rand(1)*1000);
x_i(3) =  ceil(rand(1)*1000); y_i(3) =  ceil(rand(1)*1000);
StateVector.Initial = [x_i(1),0,y_i(1),0,x_i(2),0,y_i(2),0,x_i(3),0,y_i(3),0]';
StateVector.Final = StateVector.Initial;
setappdata(handles.figure1,'StateVector',StateVector);

% --- Outputs from this function are returned to the command line.
function varargout = gui_LocalizationServer_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in StartLoclaizationButton.
function StartLoclaizationButton_Callback(hObject, eventdata, handles)
% hObject    handle to StartLoclaizationButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(handles.Flags.ComPortOpen == 0)
    try
        fopen(handles.AltSerialPort);
    catch
        set(handles.StatusText,'String','Unable to Open COM Port');
        return;
    end
    handles.Flags.ComPortOpen = 1;
end
%Reset all associated flags and states and save
handles.Flags.StartLocalization = 0;
handles.Flags.LocalizationDone = 0;
StateMachine.State.Current = 1;
StateMachine.State.Next =1;
setappdata(handles.figure1,'StateMachine',StateMachine);
guidata(hObject,handles);
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
    trailer(1) = 255 - bitand(sum(data),255); %crc check
    fwrite(handles.AltSerialPort,[header; data; trailer],'uint8');

    %Create an event each time a certain number of bytes becomes available
    %at the input buffer
    handles.Flags.StartLocalization = 1;
catch
    set(handles.StatusText,'String','Error accessing the COM port.');
    %fclose(handles.AltSerialPort);
    return;
end

%Close the COM port and set the BytesAvailable Fcn
if(handles.Flags.ComPortOpen) %ComPort open?
    fclose(handles.AltSerialPort);
end
set(handles.AltSerialPort,'BytesAvailableFcnCount',handles.DataPacketSize+6);
set(handles.AltSerialPort,'BytesAvailableFcnMode','byte')
set(handles.AltSerialPort,'BytesAvailableFcn',{@PopulateMatrix,handles});
guidata(hObject,handles);
try
    fopen(handles.AltSerialPort);
    handles.Flags.ComPortOpen = 1;
catch
    set(handles.StatusText,'String','Unable to start localization');
    return;
end
guidata(hObject,handles);
set(handles.StatusText, 'String', 'Data written successfully.');

% --- Executes on button press in StartTrackingButton.
function StartTrackingButton_Callback(hObject, eventdata, handles)
% hObject    handle to StartTrackingButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on selection change in N_nodesMenu.
function N_nodesMenu_Callback(hObject, eventdata, handles)
% hObject    handle to N_nodesMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = get(hObject,'String') returns N_nodesMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from N_nodesMenu
handles.NumberOfNodes = get(hObject,'Value') + 2;
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

BytesAvailable = get(handles.AltSerialPort,'BytesAvailable');

%Read the data to clear the buffer
if(BytesAvailable)
    A = uint8(fread(handles.AltSerialPort,handles.DataPacketSize + 6,'uint8'));
else
    return; %no data then return
end
if(handles.Flags.StartLocalization == 0) %Bogus data
    set(handles.StatusText,'String','Invalid Data');
    %     try
    %         fclose(handles.AltSerialPort);
    %     catch
    %         set(handles.StatusText,'String','Bogus data : unable to close COM port');
    %         return;
    %     end
    %     %Reenable the BytesAvailableFcn
    %     set(handles.AltSerialPort,'BytesAvailableFcnCount',handles.DataPacketSize+6);
    %     set(handles.AltSerialPort,'BytesAvailableFcnMode','byte');
    %     set(handles.AltSerialPort,'BytesAvailableFcn',{@PopulateMatrix,handles});
    %     guidata(hObject,handles);
    %     try
    %         fopen(handles.SerialPort);
    %     catch
    %         set(handles.StatusText,'String','Bogus data : unable to reopen COM port');
    %         return;
    %     end
    %     set(handles.StatusText,'String','Bogus data : Reopened COM port');
    return;
end

try
    if (A(1) ~= hex2dec('7E')) || (A(2) ~= hex2dec('01')) || ((A(3) + A(4)*256) ~= handles.DataPacketSize) || (A(4+handles.DataPacketSize+2) ~= hex2dec('FF'))
        set(handles.StatusText,'String','Error in the receieved packet form');
        fclose(handles.AltSerialPort);
        handles.Flags.ComPortOpen = 0;
        guidata(hObject,handles);
        return;
    end
    if bitand(sum(A(5:(length(A)-2))), 255) + A(length(A)-1) ~= 255
        set(handles.StatusText, 'String', 'CRC error.');
        fclose(handles.AltSerialPort);
        handles.Flags.ComPortOpen = 0;
        guidata(hObject,handles);
        return;
    end
catch
    set(handles.StatusText,'String','Error accessing the COM port');
    fclose(handles.AltSerialPort);
    handles.Flags.ComPortOpen = 0;
    guidata(hObject,handles);
    return;
end

Matrices = getappdata(handles.figure1,'LocalizationData');
for i = 0 : handles.MaxNodes - 1
    CurrentNode = i + 1;
    CurrentNodeOffset = i*handles.params_size;
    CurrentNodeStart = 5+CurrentNodeOffset;
    CurrentNodeEnd = CurrentNodeStart + handles.params_size;
    CurrentNodeData = A(CurrentNodeStart : CurrentNodeEnd);
    params = CreateParams(CurrentNodeData);
    RefNodeRadio = params.reference_node_radio;
    if RefNodeRadio == 0
        RefNodeRadio = handles.MaxNodes;
    end
    RefNodeSound = params.reference_node_sound;
    if RefNodeSound == 0
        RefNodeSound = handles.MaxNodes;
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
    Matrices.RefNode.RSSI.Radio(CurrentNode,RefNodeSound)        = params.rssi_reference_node_radio;
    Matrices.RefNode.RSSI.Sync(CurrentNode,RefNodeSound)         = params.rssi_sync;
    Matrices.RefNode.RSSI.BaseStation(CurrentNode,RefNodeSound)  = params.rssi_basestation;
    Matrices.RefNode.RSSI.Radio(RefNodeSound,CurrentNode)        = params.rssi_reference_node_radio;
    Matrices.RefNode.RSSI.Sync(RefNodeSound,CurrentNode)         = params.rssi_sync;
    Matrices.RefNode.RSSI.BaseStation(RefNodeSound,CurrentNode)  = params.rssi_basestation;
end
setappdata(handles.figure1,'LocalizationData',Matrices);
set(handles.StatusText, 'String', 'Data read successfully.');

%State dynamics
StateVector = getappdata(handles.figure1,'StateVector');
x0 = StateVector.Initial;
tspan = [0 40];
[t,x] = ode15s(@(t,x) LocalizeNewNode(t,x,handles),tspan,x0);
StateMachine  = getappdata(handles.figure1,'StateMachine');
StateMachine.State.Current = StateMachine.State.Next; %Next state
setappdata(handles.figure1,'StateMachine',StateMachine);
x_row_fin = length(x);  x_fin = x(x_row_fin,:);

%Have to save the final states (position) of each node
if(StateMachine.State.Current > 3) %Start saving only in the 3rd state
    for i = 1 : length(x_fin)/2
        XY_f(i) = x_fin(2*i - 1); %selecting only the position states
    end
    for i = 1 : length(XY_f)/2
        X_f(i) = XY_f(2*i-1);  Y_f(i) = XY_f(2*i);
    end
    if(StateMachine.State.Current == 4) %Save first three nodes
        Position.X = X_f;
        Position.Y = Y_f;
    else %Append for each node n>3
        Position = getappdata(handles.figure1,'NodePosition');
        Position.X = [Position.X, X_f]';
        Position.Y = [Position.Y, Y_f]';
    end
    %Displaying position data
    axes(handles.DisplayAxes)
    hold on,    plot(handles.DisplayAxes,Position.X,Position.Y,'o');
%     StateVector.tspan = [0 40];
%     StateVector.Final = x_fin;
%     setappdata(handles.figure1,'StateVector',StateVector);
    setappdata(handles.figure1,'NodePosition',Position);
end


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


function xdot = LocalizeNewNode(t,x,handles)
%t time span for integration
%x State Vector
%hObject    Not Passed
% eventdata  Not Passed
% handles    handles to all UI objects and UI data

%Retrieve range data
Matrices = getappdata(handles.figure1,'LocalizationData');
StateVector = getappdata(handles.figure1,'StateVector');
StateMachine = getappdata(handles.figure1,'StateMachine');
Rd = double(Matrices.Distance.Radio);
switch(StateMachine.State.Current)
    case {1,2} %Wait three states for atleast 3 nodes information to arrive
        xdot = x; %Keep initial states
        StateMachine.State.Next = StateMachine.State.Current + 1;
    case 3 % 1st 3 nodes can be localized
        if((Rd(1,2)==0)||(Rd(1,3)==0)||(Rd(2,3)==0)) %data unavailable?
            xdot = x; %Keep initial states
            StateMachine.State.Next = StateMachine.State.Current;
        else
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
            %Next state (state machine) and save initial state
            %vector for the next state (state machine)
            StateVector.Initial = [ ceil(rand(1)*1000); 0;  ceil(rand(1)*1000); 0];
            StateMachine.State.Next = StateMachine.State.Current + 1;
        end

    case 4 %4th node
        if( length(nonzeros(Rd(4,:))) < 3) % need minimum 3 distance readings
            xdot = x; %Keep initial states
            StateMachine.State.Next = StateMachine.State.Current;
        else
            %The final state vector from the previous
            %state gives the final position of the nodes
            %under consideration in that state(state machine)
            StateVector = getappdata(handles.figure1,'StateVector');
            %Retrieve position of previously localized nodes
            Position = getappdata(handles.figure1,'NodePosition');
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
            %Next state (state machine) and save initial state
            %vector for the next state (state machine)
            StateVector.Initial = [ ceil(rand(1)*1000); 0;  ceil(rand(1)*1000); 0];
            StateMachine.State.Next = StateMachine.State.Current + 1;
        end

    case {5,6,7,8} %5th 6th 7th and 8th node
        n = StateMachine.State.Current;
        if n <= handles.NumberOfNodes
            if( length(nonzeros(Rd(n,:))) < 3) % need minimum 3 distance readings
                xdot = x; %Keep initial states
                StateMachine.State.Next = StateMachine.State.Current;
            else
                %The final state vector from the previous
                %state gives the final position of the nodes
                %under consideration in that state(state machine)
                StateVector = getappdata(handles.figure1,'StateVector');
                %Retrieve position of previously localized nodes
                Position = getappdata(handles.figure1,'NodePosition');
                X_f = Position.X;   Y_f = Position.Y;
                %Renaming position variables for the nth node
                xn = x(1);      yn = x(3);
                %Renaming velocity variables
                xndot = x(2);   yndot = x(4);
                Fnx = 0; Fny = 0;
                % Gains
                k=1;            kv=1;
                for i = 1:(n-1)
                    R(n,i) = sqrt((xn-X_f(i))^2 + (yn-Y_f(i))^2);
                    Fnx = Fnx + k*(R(n,i)-Rd(n,i))*((X_f(i)-xn)/R(n,i));
                    Fny = Fny + k*(R(n,i)-Rd(n,i))*((Y_f(i)-yn)/R(n,i));
                end
                Fnx = Fnx - kv*xndot;
                Fny = Fny - kv*yndot;
                xdot=[ x(2);  % x1
                    Fnx;
                    x(4);  % y1
                    Fny];
                StateVector.Initial = [ ceil(rand(1)*1000); 0;  ceil(rand(1)*1000); 0];
                StateMachine.State.Next = StateMachine.State.Current + 1;
            end
        else
            StateMachine.State.Next = 9; %No more nodes to localize so finish
        end

    case 9 % All nodes
        %Reset the initial state vector for next localization cycle
        x_i(1) =  ceil(rand(1)*1000); y_i(1) =  ceil(rand(1)*1000);
        x_i(2) =  ceil(rand(1)*1000); y_i(2) =  ceil(rand(1)*1000);
        x_i(3) =  ceil(rand(1)*1000); y_i(3) =  ceil(rand(1)*1000);
        StateVector.Initial = [x_i(1),0,y_i(1),0,x_i(2),0,y_i(2),0,x_i(3),0,y_i(3),0]';
        StateMachine.State.Next = StateMachine.State.Current + 1;
    otherwise %Default state or state 10;
        handles.Flags.LocalizationDone = 1;
        handles.Flags.StartLocalization = 0;
end
setappdata(handles.figure1,'StateVector',StateVector);
setappdata(handles.figure1,'StateMachine',StateMachine);

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
Matrices = getappdata(handles.figure1,'LocalizationData');
BaudRate = Matrices.BaudRate;
ChosenComPort = handles.AvailableComPorts(get(handles.ComPortMenu,'Value'));
set(handles.AltSerialPort,'Port',cell2mat(ChosenComPort));
set(handles.AltSerialPort,'BaudRate',BaudRate(get(handles.BaudRateMenu,'Value')));
guidata(hObject,handles);
try
    fopen(handles.AltSerialPort);
    handles.Flags.ComPortOpen = 1;
    set(handles.StatusText,'String',['Successfully opened ',...
        get(handles.AltSerialPort,'Port'),' at ',...
        int2str(get(handles.AltSerialPort,'BaudRate')),' baud']);
catch
    set(handles.StatusText,'String',['Unable to open ',...
        get(handles.AltSerialPort,'Port'),' at ',...
        int2str(get(handles.AltSerialPort,'BaudRate')),' baud']);
    return;
end

%Save the handles structure
guidata(hObject,handles);

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
    bytes_available = get(handles.AltSerialPort,'BytesAvailable');
    if bytes_available
        A = fread(handles.AltSerialPort,bytes_available,'uint8');
        %Reading the buffer clears it
    end
    fclose(handles.AltSerialPort); %Close the port
    delete(handles.AltSerialPort); %delete serial object
    handles.AltSerialPort = 0;
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

