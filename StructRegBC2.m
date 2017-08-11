function varargout = StructRegBC2(varargin)
% STRUCTREGBC2 MATLAB code for StructRegBC2.fig
%      STRUCTREGBC2, by itself, creates a new STRUCTREGBC2 or raises the existing
%      singleton*.
%
%      H = STRUCTREGBC2 returns the handle to a new STRUCTREGBC2 or the handle to
%      the existing singleton*.
%
%      STRUCTREGBC2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STRUCTREGBC2.M with the given input arguments.
%
%      STRUCTREGBC2('Property','Value',...) creates a new STRUCTREGBC2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before StructRegBC2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to StructRegBC2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to next (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help StructRegBC2

% Last Modified by GUIDE v2.5 17-Jul-2016 19:59:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @StructRegBC2_OpeningFcn, ...
    'gui_OutputFcn',  @StructRegBC2_OutputFcn, ...
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


% --- Outputs from this function are returned to the command line.
function varargout = StructRegBC2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes just before StructRegBC2 is made visible.
function StructRegBC2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to StructRegBC2 (see VARARGIN)

% Choose default command line output for StructRegBC2
handles.output = hObject;
load('StructDef','edgesArray'); % edgesArray contains all exterior edges
load('StructRegBC1','data','DorN'); % data contains pairs of edges and their types
for i=1:length(edgesArray)
    if strcmp(data{i,2},DorN{1}) % if the edge type is Dirichlet
        handles.edgesDirichlet(i,1)=edgesArray(i); % list with the Dirichlet edges
    else
        handles.edgesNeumann(i,1)=edgesArray(i);  % list with the Neumann edges
    end
end

% Creating dataDir from scratch or reading it from StructRegBC2 file
if isfield(handles,'edgesDirichlet')==1
    handles.edgesDirichlet(handles.edgesDirichlet==0)=[];
    set(handles.listboxDir,'string',handles.edgesDirichlet);
    
    handles.dataDir=cell(length(handles.edgesDirichlet),3);
    
    if exist('./StructRegBC2.mat','file') % reuse the previous data
        load('StructRegBC2','dataDir');
        handles.dataDir=dataDir;
    else
        for i=1:length(handles.edgesDirichlet)
            handles.dataDir{i,1}=handles.edgesDirichlet(i);
            handles.dataDir{i,2}='NaN';
            handles.dataDir{i,3}='NaN';
        end
    end
    
    column1={'Boundary ID','Normal','Tangential'};
    uitable('units','Normalized','Position',[0.475,0.10,0.13,0.44],...
        'Data',handles.dataDir,...
        'ColumnName',column1,'RowName',[]);
end

% Creating dataNeu from scratch or reading it from StructRegBC2 file
if isfield(handles,'edgesNeumann')==1

    handles.edgesNeumann(handles.edgesNeumann==0)=[];
    set(handles.listboxNeu,'string',handles.edgesNeumann);
    
    handles.dataNeu=cell(length(handles.edgesNeumann),3);
    
    if exist('./StructRegBC2.mat','file') % reuse the previous data
        load('StructRegBC2','dataNeu');
        handles.dataNeu=dataNeu;
    else
        for i=1:length(handles.edgesNeumann)
            handles.dataNeu{i,1}=handles.edgesNeumann(i);
            handles.dataNeu{i,2}='NaN';
            handles.dataNeu{i,3}='NaN';
        end
    end
    
    column2={'Boundary ID','Normal','Tangential'};
    uitable('units','Normalized','Position',[0.845,0.10,0.13,0.44],...
        'Data',handles.dataNeu,...
        'ColumnName',column2,'RowName',[]);
end

% ****** DRAW THE MESH ******

%load the array of edges
load('StructDef','L','B','Nx','Ny');
%load the mesh data
load('StructRegBC1','nodes','edges_nodes','edges_loops','loops_nodes',...
    'loops_edges');

nel = Nx*Ny ;        % Total Number of Elements in the Mesh
nnel = 4 ;           % Number of nodes per Element

% Number of points on the Length and Width
npx = Nx+1 ;
npy = Ny+1 ;
nnode = npx*npy ;      % Total Number of Nodes in the Mesh

% Number of edges on the Length and Width
nedgex = Nx*npy;       % Number of Horizontal Edges
nedgey = Ny*npx;       % Number of Vertical Edges
nedge = nedgex+nedgey; % Total Number of Edges in the Mesh

% For drawing purposes
limxmin = min(nodes(:,1));
limxmax = max(nodes(:,1));
limymin =  min(nodes(:,2));
limymax =  max(nodes(:,2));

%
% Plotting the Finite Element Mesh
% Initialization of the required matrices
X = zeros(nnel,nel) ;
Y = zeros(nnel,nel) ;
% Extract X,Y coordinates for the (iel)-th element
for iel = 1:nel
    X(:,iel) = nodes(loops_nodes(iel,:),1) ;
    Y(:,iel) = nodes(loops_nodes(iel,:),2) ;
end

patch(X,Y,'w');
axis([limxmin-0.01*abs(limxmin) limxmax+0.01*abs(limxmax) limymin-0.01*abs(limymin) limymax+0.01*abs(limymax)]);
axis equal;
axis off ;

% To display Node Numbers % Element Numbers
axpos = getpixelposition(handles.axes3); % position & dimension of the axes object
% Define button's weight and height
bweight = 60;
bheight = 20;
pos = [((axpos(1)+axpos(3))/2) (axpos(2)-1.5*bheight) bweight bheight]; % align the second button with the center of the axes obj limit
ShowNodes = uicontrol('style','toggle','string','Nodes',....
    'position',[(pos(1)-2*bweight) pos(2) pos(3) pos(4)],'background',...
    'white','units','Normalized');

ShowEdges = uicontrol('style','toggle','string','Edges',....
    'position',[pos(1) pos(2) pos(3) pos(4)],'background','white',...
    'units','Normalized');

ShowElements = uicontrol('style','toggle','string','Elements',....
    'position',[(pos(1)+2*bweight) pos(2) pos(3) pos(4)],'background',...
    'white','units','Normalized');

set(ShowNodes,'callback',...
    {@SHOWNODES,ShowEdges,ShowElements,nodes,edges_nodes,loops_nodes,X,Y,nnode,nedge,nel});
set(ShowEdges,'callback',...
    {@SHOWEDGES,ShowNodes,ShowElements,nodes,edges_nodes,loops_nodes,X,Y,nnode,nedge,nel});
set(ShowElements,'callback',....
    {@SHOWELEMENTS,ShowNodes,ShowEdges,nodes,edges_nodes,loops_nodes,X,Y,nnode,nedge,nel});

%-----------------------------------------------------------------

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in next.
function next_Callback(hObject, eventdata, handles)
% hObject    handle to next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

NaND = 0; NaNN = 0;

if isfield(handles,'edgesDirichlet')
    edgesDirichlet=handles.edgesDirichlet;
    dataDir=handles.dataDir;
    % Converting strings in data arrays to matrices
    D = cellfun(@str2num,dataDir(:,2:end),'UniformOutput',0);
    % Looking for NaN in D & N
    NaND = cellfun(@any,cellfun(@isnan,D,'UniformOutput',0));
end
if isfield(handles,'edgesNeumann')
    edgesNeumann=handles.edgesNeumann;
    dataNeu=handles.dataNeu;
    % Converting strings in data arrays to matrices
    N = cellfun(@str2num,dataNeu(:,2:end),'UniformOutput',0);
    % Looking for NaN in D & N
    NaNN = cellfun(@any,cellfun(@isnan,N,'UniformOutput',0));
end

%%
% Checking the data

if any(all(NaND,2)) % if all Dirichlet conditions on an edge are NaN
    errordlg('All Dirichlet boundary conditions are set to NaN at least for one boundary.','Invalid input','modal');
    return;
end

if any(any(NaNN,2)) % if any Neumann conditions are NaN
    errordlg('At least one Neumann boundary condition is set to NaN.','Invalid input','modal');
    return;
end

%%
% saving the workspace
if ~exist('./StructRegBC2.mat','file')
    dummy = 0;
    save('StructRegBC2','dummy'); % if save file does not exist, it creates it
end

load('StructDef','DirName','FileName');
% if requested by the user, saves to file

if isfield(handles,'edgesDirichlet')
    if ~isempty(DirName)
        save(fullfile(DirName,FileName),'-append',...
            'edgesDirichlet','dataDir');
    end
    save('StructRegBC2','-append','edgesDirichlet','dataDir');
end
if isfield(handles,'edgesNeumann')
    if ~isempty(DirName)
        save(fullfile(DirName,FileName),'-append',...
            'edgesNeumann','dataNeu');
    end
    save('StructRegBC2','-append','edgesNeumann','dataNeu');
end

close(handles.figure1); % closing the GUI window

% Trying to close Visualize if it was opened
try
    close('Visualize');
catch
end

CheckStructReg;


% --- Executes on button press in assignForce.
function assignForce_Callback(hObject, eventdata, handles)
% hObject    handle to assignForce (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

itemlist = get(handles.listboxNeu,'Value'); % list of selected items in the listbox
nitems = length(itemlist); % number of selected items in the listbox
ForceDirection = get(handles.popupmenu_NTNeu,'Value');

for ii = 1:nitems
    crtitem = itemlist(ii);
    handles.dataNeu{crtitem,ForceDirection+1}=get(handles.editForce,'string');
end

guidata(hObject,handles);

column2={'Boundary ID','Normal','Tangential'};
uitable('units','Normalized','Position',[0.845,0.10,0.13,0.44],...
    'Data',handles.dataNeu,...
    'ColumnName',column2,'RowName',[]);


% --- Executes on button press in assignDisp.
function assignDisp_Callback(hObject, eventdata, handles)
% hObject    handle to assignDisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

itemlist = get(handles.listboxDir,'Value'); % list of selected items in the listbox
nitems = length(itemlist); % number of selected items in the listbox
DispDirection = get(handles.popupmenu_NTDir,'Value');

for ii = 1:nitems
    crtitem = itemlist(ii);
    handles.dataDir{crtitem,DispDirection+1}=get(handles.editDisp,'string');
end

guidata(hObject,handles);

column1={'Boundary ID','Normal','Tangential'};
uitable('units','Normalized','Position',[0.475,0.10,0.13,0.44],...
    'Data',handles.dataDir,...
    'ColumnName',column1,'RowName',[]);


% --- Executes on button press in resetNeumann.
function resetNeumann_Callback(hObject, eventdata, handles)
% hObject    handle to resetNeumann (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
for i=1:length(handles.edgesNeumann)
    handles.dataNeu{i,2}='NaN';
    handles.dataNeu{i,3}='NaN';
end
guidata(hObject,handles);

column2={'Boundary ID','Normal','Tangential'};
uitable('units','Normalized','Position',[0.845,0.10,0.13,0.44],...
    'Data',handles.dataNeu,...
    'ColumnName',column2,'RowName',[]);



% --- Executes on button press in resetDirichlet.
function resetDirichlet_Callback(hObject, eventdata, handles)
% hObject    handle to resetDirichlet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
for i=1:length(handles.edgesDirichlet)
    handles.dataDir{i,2}='NaN';
    handles.dataDir{i,3}='NaN';
end
guidata(hObject,handles);

column1={'Boundary ID','Normal','Tangential'};
uitable('units','Normalized','Position',[0.475,0.10,0.13,0.44],...
    'Data',handles.dataDir,...
    'ColumnName',column1,'RowName',[]);


% --- Executes on button press in previous.
function previous_Callback(hObject, eventdata, handles)
% hObject    handle to previous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(StructRegBC2);

% Trying to close Visualize if it was opened
try
    close('Visualize');
catch
end

StructRegBC1;


% --- Executes on selection change in listboxNeu.
function listboxNeu_Callback(hObject, eventdata, handles)
% hObject    handle to listboxNeu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxNeu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxNeu


% --- Executes during object creation, after setting all properties.
function listboxNeu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxNeu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editForce_Callback(hObject, eventdata, handles)
% hObject    handle to editForce (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editForce as text
%        str2double(get(hObject,'String')) returns contents of editForce as a double

[Flux, status] = str2num(get(hObject,'string'));
if any(isnan(Flux)) || ~status  % if the input is something else than
                                     % a vector of reals
    set(hObject,'String','');
    errordlg('Flux field must have real value','Invalid input','modal');
    uicontrol(hObject);
    return;
end


% --- Executes during object creation, after setting all properties.
function editForce_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editForce (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listboxDir.
function listboxDir_Callback(hObject, eventdata, handles)
% hObject    handle to listboxDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxDir contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxDir


% --- Executes during object creation, after setting all properties.
function listboxDir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editDisp_Callback(hObject, eventdata, handles)
% hObject    handle to editDisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDisp as text
%        str2double(get(hObject,'String')) returns contents of editDisp as a double
[U, status] = str2num(get(hObject,'string'));
if ~status  % if the input is something else than a vector of reals (or a NaN)
    set(hObject,'String','');
    errordlg('Sate field must have real value','Invalid input','modal');
    uicontrol(hObject);
    return;
end

if length(U)~=1 && any(isnan(U)) % if more than one term in input and one of them is NaN
    set(hObject,'string','');
    errordlg('A NaN Dirichlet boundary condition cannot be defined along with any other term','Invalid input','modal');
    uicontrol(hObject);
    return;
end

% --- Executes during object creation, after setting all properties.
function editDisp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_NTDir.
function popupmenu_NTDir_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_NTDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_NTDir contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_NTDir


% --- Executes during object creation, after setting all properties.
function popupmenu_NTDir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_NTDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_NTNeu.
function popupmenu_NTNeu_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_NTNeu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_NTNeu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_NTNeu


% --- Executes during object creation, after setting all properties.
function popupmenu_NTNeu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_NTNeu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in enlarge.
function enlarge_Callback(hObject, eventdata, handles)
% hObject    handle to enlarge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

VisualizeReg;
