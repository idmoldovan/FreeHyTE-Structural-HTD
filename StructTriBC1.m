%% STANDARD, PREDEFINED FUNCTIONS
% No changes were made to the following functions, predefined by GUIDE

function varargout = StructTriBC1(varargin)
% STRUCTTRIBC1 MATLAB code for StructTriBC1.fig
%      STRUCTTRIBC1, by itself, creates a new STRUCTTRIBC1 or raises the existing
%      singleton*.
%
%      H = STRUCTTRIBC1 returns the handle to a new STRUCTTRIBC1 or the handle to
%      the existing singleton*.
%
%      STRUCTTRIBC1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STRUCTTRIBC1.M with the given input arguments.
%
%      STRUCTTRIBC1('Property','Value',...) creates a new STRUCTTRIBC1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before StructTriBC1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to StructTriBC1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help StructTriBC1

% Last Modified by GUIDE v2.5 03-Jun-2017 13:50:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @StructTriBC1_OpeningFcn, ...
                   'gui_OutputFcn',  @StructTriBC1_OutputFcn, ...
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
function varargout = StructTriBC1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%% OPENING FUNCTION
% * Executes just before StructTriBC1 is made visible;
% * Reads the mesh data from the |mat| file of the previous GUI and
% constructs the topological matrices;
% * If the |mat| file associated to the current GUI exists (meaning that
% the previous GUI was not changed), it loads the boundary information,
% otherwise it sets all boundaries to Dirichlet.
function StructTriBC1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to StructTriBC1 (see VARARGIN)

%%
% Creates the |handles| structure
% Choose default command line output for |StructTriBC1|
handles.output = hObject;

%%
% Loading the mesh information (is different for Regular and Triangular
% routines, but essentially does the same thing).
load('StructDef','p','t');

%%
% Calls an external function that creates all topological matrices. In the
% Regular routine, this is done in the opening function, rather than
% externally.
[nodes, edges_nodes, edges_loops, loops_nodes, loops_edges] = CreateEdgeLoop(p,t);

%%
% Getting the exterior edges. In the Regular routine, this information is
% obtained from the |mat| file of the previous GUI.
edgesArray = find(edges_loops(:,2)==0);

%%
% Publishing the topological information in the |handles| structure, to
% make it available to all functions associated to the current GUI.
set(handles.listbox1,'string',edgesArray);
handles.DorN = get(handles.popupmenu1,'String');
handles.data = cell(length(edgesArray),2);
handles.edgesArray = edgesArray;
handles.nodes = nodes;
handles.edges_nodes = edges_nodes;
handles.edges_loops = edges_loops;
handles.loops_nodes = loops_nodes;
handles.loops_edges = loops_edges;

%%
% If there exists a local |mat| file, it loads its key elements to fill in
% the fields with data taken from the previous run...
if exist('./StructTriBC1.mat','file') % reuse the previous data
    load('StructTriBC1','data');
    handles.data=data;
%%
% ... otherwise it just presets all boundaries to 'Dirichlet' type.
else
    for i=1:length(edgesArray)
        handles.data{i,1}=edgesArray(i);
        handles.data{i,2}=handles.DorN{1};  % Dirichlet as predefined boundary type
    end
end

%%
% Creates the table where the boundaries are listed, along with their types
column = {'Boundary ID','Boundary Type'};
uitable('units','Normalized','Position',[0.76, 0.17, 0.13, 0.35],'Data',...
    handles.data,'ColumnName',column,'RowName',[]);

%%
% Generates the code for drawing the mesh, along with the mesh information
% buttons

nel = size(loops_nodes,1);        % Total Number of Elements in the Mesh
nnel = size(loops_nodes,2);           % Number of nodes per Element
nnode = size(nodes,1);      % Total Number of Nodes in the Mesh
nedge = size(edges_nodes,1); % Total Number of Edges in the Mesh

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

%%
% Updates the handles structure
guidata(hObject, handles);


%% NEXT BUTTON
% * Executes on button press in |pushbutton_next|;
% * It recovers all data provided by the user in the GUI fields and the
% topological matrices computed in the opening function;
% * It checks for relevant (i.e. mesh and boundary type) changes as 
% compared to the previous |mat| file. If such changes exist, it deletes 
% the |mat| file corresponding to the next GUI to force its definition from 
% scratch. If no mesh changes were detected, the next GUI will load its 
% previous version;
% * If the user asked for the model to be saved, it saves the information
% regarding this GUI in the specified file and folder.
function pushbutton_next_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
% Recovering the data created or set by the user in the current GUI
data=handles.data;
DorN=handles.DorN;
edgesArray = handles.edgesArray;
nodes = handles.nodes;
edges_nodes = handles.edges_nodes;
edges_loops = handles.edges_loops;
loops_nodes = handles.loops_nodes;
loops_edges = handles.loops_edges;

%%
% Check if critical changes were made in the current GUI session. 
if exist('./StructTriBC1.mat','file')
    save('ToDetectChanges','data'); % temporary file to detect critical changes
    % Fields whose change triggers the change of the mesh
    CriticalFields = {'data'};
    % Checks the old StructRegBC1 and the new ToDetectChanges for changes in
    % the critical fields. ChangesQ = 1 if there are changes, 0 otherwise
    ChangesQ = any(~isfield(comp_struct(load('StructTriBC1.mat'),...
        load('ToDetectChanges.mat')),CriticalFields));
    % deleting the auxiliary file
    delete('ToDetectChanges.mat');
%%
% If the |mat| file associated to the current GUI session does not exist,
% it is likely because the previous GUI was modified. In this case the
% next GUI needs to be reinitialized.
else
    ChangesQ = 1;
end

%%
% Saving the workspace to the local |mat| file. If requested by the
% user, the GUI information is also _appended_ to the file specified by the
% user.
load('StructDef','DirName','FileName');
% if requested by the user, saves to file
if ~isempty(DirName)
    save(fullfile(DirName,FileName),'-append',...
        'data','DorN','nodes','edges_nodes','edges_loops',...
        'loops_nodes','loops_edges','edgesArray');
end

% saves the local data file
save('StructTriBC1',...
    'data','DorN','nodes','edges_nodes','edges_loops',...
    'loops_nodes','loops_edges','edgesArray');

% saves the data required by VisualizeTri to the local data file
save('Visualize',...
    'nodes','edges_nodes','edges_loops',...
    'loops_nodes','loops_edges');

%%
% Preparing to exit.
close(handles.figure1);

% Trying to close Visualize if it was opened
try
    close('Visualize');
catch
end

%%
% If there are relevant changes, it physically deletes the |mat| files
% associated to the next GUI.
if ChangesQ
    delete('StructTriBC2.mat');
end
StructTriBC2;


%% RESET BUTTON
% * Executes on button press in |pushbutton_reset|.
% * It resets all tables and defines all boundaries as Dirichlet.
function pushbutton_reset_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
% Setting all edges to Dirichlet type...
edgesArray = handles.edgesArray;
handles.DorN = get(handles.popupmenu1,'String');
handles.data = cell(length(edgesArray),2);
for i=1:length(edgesArray)
    handles.data{i,1}=edgesArray(i); 
    handles.data{i,2}=handles.DorN{1};  % Dirichlet as predefined boundary type
end

%%
% ... and redrawing the table.
column = {'Boundary ID','Boundary Type'};
uitable('units','Normalized','Position',[0.76, 0.17, 0.13, 0.35],'Data',...
    handles.data,'ColumnName',column,'RowName',[]);


%% PREVIOUS BUTTON
% * Executes on button press in |pushbutton_previous|;
% * Just closes the current GUI and launches the previous one. All changes
% made in the current GUI are lost.
function pushbutton_previous_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_previous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(StructTriBC1);

% Trying to close Visualize if it was opened
try
    close('Visualize');
catch
end

StructDef(2);



%% ASSIGN BUTTON
% * Executes on button press in |assign|;
% * Fills in the boundary type table for the boundaries selected in
% |listbox1| with 'Dirichlet' or 'Neumann', as defined in |popupmenu1|.
function assign_Callback(hObject, eventdata, handles)
% hObject    handle to assign (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
% Gets the selected items in |listbox1|
itemlist = get(handles.listbox1,'value'); % list of selected items in the listbox
nitems = length(itemlist); % number of selected items in the listbox

%%
% Sets the |data| for the selected items to the type selected in
% |popupmenu1|.
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

%%
% Updates the handles structure
guidata(hObject,handles);

%% 
% Redraws the table
column = {'Boundary ID','Boundary Type'};
uitable('units','Normalized','Position',[0.76, 0.17, 0.13, 0.35],'Data',...
    handles.data,'ColumnName',column,'RowName',[]);


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

% saves the data required by Visualize to the local data file
save('Visualize',...
    'nodes','edges_nodes','edges_loops',...
    'loops_nodes','loops_edges');

VisualizeTri;



%% GUI ENTITIES GENERATION CODE
% * Automatically generated code
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

