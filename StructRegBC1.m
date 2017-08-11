function varargout = StructRegBC1(varargin)
% STRUCTREGBC1 MATLAB code for StructRegBC1.fig
%      STRUCTREGBC1, by itself, creates a new STRUCTREGBC1 or raises the existing
%      singleton*.
%
%      H = STRUCTREGBC1 returns the handle to a new STRUCTREGBC1 or the handle to
%      the existing singleton*.
%
%      STRUCTREGBC1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STRUCTREGBC1.M with the given input arguments.
%
%      STRUCTREGBC1('Property','Value',...) creates a new STRUCTREGBC1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before StructRegBC1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to StructRegBC1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help StructRegBC1

% Last Modified by GUIDE v2.5 15-Jul-2016 10:27:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @StructRegBC1_OpeningFcn, ...
                   'gui_OutputFcn',  @StructRegBC1_OutputFcn, ...
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
function varargout = StructRegBC1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes just before StructRegBC1 is made visible.
function StructRegBC1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to StructRegBC1 (see VARARGIN)

% Choose default command line output for StructRegBC1
handles.output = hObject;
%load the array of edges
load('StructDef','edgesArray','L','B','Nx','Ny');

set(handles.listbox1,'string',edgesArray);
handles.DorN = get(handles.popupmenu1,'String');
handles.data = cell(length(edgesArray),2);

if exist('./StructRegBC1.mat','file') % reuse the previous data
    load('StructRegBC1','data');
    handles.data=data;
else
    for i=1:length(edgesArray)
        handles.data{i,1}=edgesArray(i);
        handles.data{i,2}=handles.DorN{1};  % Dirichlet as predefined boundary type
    end 
end


column = {'Boundary ID','Boundary Type'};
uitable('units','Normalized','Position',[0.76, 0.17, 0.13, 0.35],'Data',...
    handles.data,'ColumnName',column,'RowName',[]);


%----Generate Mesh Code-----------------------------------------------

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

% Discretizing the Length and Width of the plate
nx = linspace(0,L,npx) ;
ny = linspace(0,B,npy) ;
[xx,yy] = meshgrid(nx,ny) ;
% To get the Nodal Connectivity Matrix
nodes = [xx(:) yy(:)] ;
NodeMap = 1:nnode ;
loops_nodes = zeros(nel,nnel) ;
% If elements along the X-axes and Y-axes are equal
if npx==npy
    NodeMap = reshape(NodeMap,npx,npy);
    loops_nodes(:,1) = reshape(NodeMap(1:npx-1,1:npy-1),nel,1);
    loops_nodes(:,2) = reshape(NodeMap(2:npx,1:npy-1),nel,1);
    loops_nodes(:,3) = reshape(NodeMap(2:npx,2:npy),nel,1);
    loops_nodes(:,4) = reshape(NodeMap(1:npx-1,2:npy),nel,1);
% If the elements along the axes are different
else%if npx>npy
    NodeMap = reshape(NodeMap,npy,npx);
    loops_nodes(:,1) = reshape(NodeMap(1:npy-1,1:npx-1),nel,1);
    loops_nodes(:,2) = reshape(NodeMap(2:npy,1:npx-1),nel,1);
    loops_nodes(:,3) = reshape(NodeMap(2:npy,2:npx),nel,1);
    loops_nodes(:,4) = reshape(NodeMap(1:npy-1,2:npx),nel,1);
end

% Obtaining the edge matrix

% Finding the nodes for each edge
edgesx=[reshape(NodeMap(1:npy-1,:),nedgey,1),...
    reshape(NodeMap(2:npy,:),nedgey,1)];
edgesx(1:Ny,:)=fliplr(edgesx(1:Ny,:));
edgesy=[reshape(NodeMap(:,2:npx),nedgex,1),...
    reshape(NodeMap(:,1:npx-1),nedgex,1)];
edgesy(1:npy:end,:)=fliplr(edgesy(1:npy:end,:));
edges_nodes = [edgesx ; edgesy];
%Finding the neighbouring elements for each edge
edges_loops = zeros(nedge,2);
for i=1:nedge
    [row1, column1] = find(loops_nodes==edges_nodes(i,1));
    [row2, column2] = find(loops_nodes==edges_nodes(i,2));
    edges_loops(i,:)=intersect(row1,row2)';
    if edges_loops(i,1)==edges_loops(i,2) % a single neighbour
        edges_loops(i,2)=0;              % right neighbour is set to zero
    elseif edges_loops(i,1) > edges_loops(i,2) % left neighbour larger than right
        edges_loops(i,:)=fliplr(edges_loops(i,:)); % flip them (left always
    end                                 % smallest in regular meshes)
end
%Finding the edges of each finite element
loops_edges = zeros(nel,4);
for i=1:nel
    [row1, column1] = find(edges_loops==i);
    loops_edges(i,:)=row1';
end

% storing all structures
handles.edgesArray = edgesArray;
handles.nodes = nodes;
handles.edges_nodes = edges_nodes;
handles.edges_loops = edges_loops;
handles.loops_nodes = loops_nodes;
handles.loops_edges = loops_edges;


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

patch(X,Y,'w')
axis([limxmin-0.01*abs(limxmin) limxmax+0.01*abs(limxmax) limymin-0.01*abs(limymin) limymax+0.01*abs(limymax)]);
axis equal;
axis off ;

% To display Node Numbers % Element Numbers
axpos = getpixelposition(handles.axes2); % position & dimension of the axes object
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
%setFigDockGroup(fh,'handles.axes1')
% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in pushbutton_next.
function pushbutton_next_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=handles.data;
DorN=handles.DorN;
edgesArray = handles.edgesArray;
nodes = handles.nodes;
edges_nodes = handles.edges_nodes;
edges_loops = handles.edges_loops;
loops_nodes = handles.loops_nodes;
loops_edges = handles.loops_edges;

% check if critical changes were made in the current GUI session
if exist('./StructRegBC1.mat','file')
    save('ToDetectChanges','data'); % temporary file to detect critical changes
    % Fields whose change triggers the change of the mesh
    CriticalFields = {'data'};
    % Checks the old StructRegBC1 and the new ToDetectChanges for changes in
    % the critical fields. ChangesQ = 1 if there are changes, 0 otherwise
    ChangesQ = any(~isfield(comp_struct(load('StructRegBC1.mat'),...
        load('ToDetectChanges.mat')),CriticalFields));
    % deleting the auxiliary file
    delete('ToDetectChanges.mat');
else
    ChangesQ = 1;
end
    
% saving the workspace
load('StructDef','DirName','FileName');
% if requested by the user, saves to file
if ~isempty(DirName)
    save(fullfile(DirName,FileName),'-append',...
        'data','DorN','nodes','edges_nodes','edges_loops',...
        'loops_nodes','loops_edges');
end

% saves the local data file
save('StructRegBC1','data','DorN','nodes','edges_nodes','edges_loops',...
    'loops_nodes','loops_edges');

% saves the data required by VisualizeReg to the local data file
save('Visualize',...
    'nodes','edges_nodes','edges_loops',...
    'loops_nodes','loops_edges');

close(handles.figure1);

% Trying to close Visualize if it was opened
try
    close('Visualize');
catch
end

if ChangesQ          % if changes were detected
    delete('StructRegBC2.mat');
end
StructRegBC2; 


% --- Executes on button press in pushbutton_reset.
function pushbutton_reset_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load('StructDef','edgesArray');
handles.DorN = get(handles.popupmenu1,'String');
handles.data = cell(length(edgesArray),2);
for i=1:length(edgesArray)
    handles.data{i,1}=edgesArray(i); 
    handles.data{i,2}=handles.DorN{1};  % Dirichlet as predefined boundary type
end

column = {'Boundary ID','Boundary Type'};
uitable('units','Normalized','Position',[0.76, 0.17, 0.13, 0.35],'Data',...
    handles.data,'ColumnName',column,'RowName',[]);



% --- Executes on button press in pushbutton_previous.
function pushbutton_previous_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_previous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(StructRegBC1);

% Trying to close Visualize if it was opened
try
    close('Visualize');
catch
end

StructDef(2);



% --- Executes on button press in assign.
function assign_Callback(hObject, eventdata, handles)
% hObject    handle to assign (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

itemlist = get(handles.listbox1,'value'); % list of selected items in the listbox
nitems = length(itemlist); % number of selected items in the listbox

if get(handles.popupmenu1,'value')==1   % if Dirichlet is selected
    for ii = 1:nitems
        crtitem = itemlist(ii);
        handles.data{crtitem,2}=handles.DorN{1};
    end
else                                    % if Neumann is selected
    for ii = 1:nitems
        crtitem = itemlist(ii);
        handles.data{crtitem,2}=handles.DorN{2};
    end
end
guidata(hObject,handles);

column = {'Boundary ID','Boundary Type'};
uitable('units','Normalized','Position',[0.76, 0.17, 0.13, 0.35],'Data',...
    handles.data,'ColumnName',column,'RowName',[]);



% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


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


% --- Executes on button press in enlarge.
function enlarge_Callback(hObject, eventdata, handles)
% hObject    handle to enlarge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%%
% Recovering the data created or set by the user in the current GUI
nodes = handles.nodes;
edges_nodes = handles.edges_nodes;
edges_loops = handles.edges_loops;
loops_nodes = handles.loops_nodes;
loops_edges = handles.loops_edges;

% saves the data required by VisualizeReg to the local data file
save('Visualize',...
    'nodes','edges_nodes','edges_loops',...
    'loops_nodes','loops_edges');

VisualizeReg;
