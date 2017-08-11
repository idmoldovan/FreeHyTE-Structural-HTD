%% STANDARD, PREDEFINED FUNCTIONS
% No changes were made to the following functions, predefined by GUIDE

function varargout = StructTriBC2(varargin)
% STRUCTTRIBC2 MATLAB code for StructTriBC2.fig
%      STRUCTTRIBC2, by itself, creates a new STRUCTTRIBC2 or raises the existing
%      singleton*.
%
%      H = STRUCTTRIBC2 returns the handle to a new STRUCTTRIBC2 or the handle to
%      the existing singleton*.
%
%      STRUCTTRIBC2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STRUCTTRIBC2.M with the given input arguments.
%
%      STRUCTTRIBC2('Property','Value',...) creates a new STRUCTTRIBC2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before StructTriBC2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to StructTriBC2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to next (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help StructTriBC2

% Last Modified by GUIDE v2.5 15-Jul-2016 10:27:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @StructTriBC2_OpeningFcn, ...
    'gui_OutputFcn',  @StructTriBC2_OutputFcn, ...
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
function varargout = StructTriBC2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%% OPENING FUNCTION
% * Executes just before |StructTriBC2| is made visible;
% * Reads the mesh data from the |mat| file of the previous GUI and
% constructs the lists of Dirichlet and Neumann boundaries;
% * If the |mat| file of the current GUI exists (meaning that the 
% previous GUI was not changed), it loads the boundary information,
% otherwise it sets all boundaries to Dirichlet.

function StructTriBC2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to StructTriBC2 (see VARARGIN)

%%
% Creates the |handles| structure
% Choose default command line output for |StructTriBC2|
handles.output = hObject;

%%
% Loads the exterior edges and their types from the previous GUI
load('StructTriBC1','edgesArray','data','DorN'); % data contains pairs of edges and their types

%%
% Creates two arrays listing the exterior Dirichlet and Neumann edges
% (|edgesDirichlet| and |edgesNeumann|)
for i=1:length(edgesArray)
    if strcmp(data{i,2},DorN{1}) % if the edge type is Dirichlet
        handles.edgesDirichlet(i,1)=edgesArray(i); % list with the Dirichlet edges
    else
        handles.edgesNeumann(i,1)=edgesArray(i);  % list with the Neumann edges
    end
end

%%
% *Dirichlet boundary conditions*

%%
% Generates the |dataDir| matrix to store the id of the exterior Dirichlet 
% edges and their enforced boundary conditions in the normal and tangential
% directions. The |dataDir| matrix is created from from scratch or imported
% from the |StructRegBC2| file, if such file exists (meaning that the
% previous GUI was left unchanged).
% It is recalled that if a displacement boundary condition is enforced in a
% single direction and left free in the other, the boundary should be
% defined as Dirichlet and the kinematic boundary condition in the free
% direction should be defined as NaN.

if isfield(handles,'edgesDirichlet')==1
    handles.edgesDirichlet(handles.edgesDirichlet==0)=[];
    
    %%
    % Writes all |edgesDirichlet| into the |listboxDir|
    set(handles.listboxDir,'string',handles.edgesDirichlet);
    
    %%
    % Creates the |dataDir| matrix.
    handles.dataDir=cell(length(handles.edgesDirichlet),3);
    
    %%
    % If there exists a local |mat| file, it loads it to fill in
    % the fields with data taken from the previous next...
    if exist('./StructTriBC2.mat','file') % reuse the previous data
        load('StructTriBC2','dataDir');
        handles.dataDir=dataDir;
    %%
    % ... otherwise it just sets all Dirichlet boundary conditions to NaN    
    else
        for i=1:length(handles.edgesDirichlet)
            handles.dataDir{i,1}=handles.edgesDirichlet(i);
            handles.dataDir{i,2}='NaN';
            handles.dataDir{i,3}='NaN';
        end
    end
    
    %%
    column1={'Boundary ID','Normal','Tangential'};
    uitable('units','Normalized','Position',[0.475,0.10,0.13,0.44],...
        'Data',handles.dataDir,...
        'ColumnName',column1,'RowName',[]);
end

%%
% *Neumann boundary conditions*

%%
% Generates the |dataNeu| matrix to store the id of the exterior Neumann 
% edges and their enforced boundary conditions in the normal and tangential
% directions. The |dataNeu| matrix is created from from scratch or imported
% from the |StructRegBC2| file, if such file exists (meaning that the
% previous GUI was left unchanged).

if isfield(handles,'edgesNeumann')==1

    handles.edgesNeumann(handles.edgesNeumann==0)=[];
    
    %%
    % Writes all |edgesNeumann| into the |listboxNeu|
    set(handles.listboxNeu,'string',handles.edgesNeumann);
    
    %%
    % Creates the |dataNeu| matrix.
    handles.dataNeu=cell(length(handles.edgesNeumann),3);
    
    
    %%
    % If there exists a local |mat| file, it loads it to fill in
    % the fields with data taken from the previous next...
    if exist('./StructTriBC2.mat','file')  % reuse the previous data
        load('StructTriBC2','dataNeu');
        handles.dataNeu=dataNeu;
    %%
    % ... otherwise it just sets all Neumann boundary conditions to NaN
    else
        for i=1:length(handles.edgesNeumann)
            handles.dataNeu{i,1}=handles.edgesNeumann(i);
            handles.dataNeu{i,2}='NaN';
            handles.dataNeu{i,3}='NaN';
        end
    end
    
    %%
    % Creates the table where the Neumann boundary conditions are listed, 
    % along with their types
    column2={'Boundary ID','Normal','Tangential'};
    uitable('units','Normalized','Position',[0.845,0.10,0.13,0.44],...
        'Data',handles.dataNeu,...
        'ColumnName',column2,'RowName',[]);
end


%%
% Generates the code for drawing the mesh, along with the mesh information
% buttons

%load the mesh data
load('StructTriBC1','nodes','edges_nodes','edges_loops','loops_nodes',...
    'loops_edges');

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

patch(X,Y,'w')
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


%%
% Updates the handles structure
guidata(hObject, handles);

%% NEXT FUNCTION
% * Executes on button press in |next|;
% * Reads the boundary condition data, stores it in the local |mat| file
% and starts the checking GUI;
% * If the user asked for the model to be saved, it saves the information
% regarding this GUI in the specified file and folder.

function next_Callback(hObject, eventdata, handles)
% hObject    handle to next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
% Recovering the user-defined data 
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
% Saving the local data to save files
% If StructTriBC2.mat does not exist, it creates it
if ~exist('./StructTriBC2.mat','file')
    dummy = 0;
    save('StructTriBC2','dummy');
end

%%
% If requested by the user, appends the local data to the save file
load('StructDef','DirName','FileName');

if isfield(handles,'edgesDirichlet')
    if ~isempty(DirName)
        save(fullfile(DirName,FileName),'-append',...
            'edgesDirichlet','dataDir');
    end
    save('StructTriBC2','-append','edgesDirichlet','dataDir');
end
if isfield(handles,'edgesNeumann')
    if ~isempty(DirName)
        save(fullfile(DirName,FileName),'-append',...
            'edgesNeumann','dataNeu');
    end
    save('StructTriBC2','-append','edgesNeumann','dataNeu');
end  

%%
% Closes everything and launches the ckecking GUI
close(handles.figure1); % closing the GUI window

% Trying to close Visualize if it was opened
try
    close('Visualize');
catch
end

CheckStructTri;


%% ASSIGN FORCE FUNCTION
% * Executes on button press in |assignForce|;
% * Fills in the enforced load table for the boundaries selected in
% |listboxNeu| with the string defined in |editForce|, in the direction
% defined in the |popupmenu_NTNeu|.

function assignForce_Callback(hObject, eventdata, handles)
% hObject    handle to assignForce (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
% Gets the selected boundaries in |listboxNeu| and the direction of the
% load in |popupmenu_NTNeu|
itemlist = get(handles.listboxNeu,'Value'); % list of selected items in the listbox
nitems = length(itemlist); % number of selected items in the listbox
ForceDirection = get(handles.popupmenu_NTNeu,'Value');

%%
% Writes the force definition in the column of |dataNeu| corresponding to
% the direction the force is applied into.
for ii = 1:nitems
    crtitem = itemlist(ii);
    handles.dataNeu{crtitem,ForceDirection+1}=get(handles.editForce,'string');
end

%%
% Updates the handles structure
guidata(hObject,handles);

%%
% Redraws the table where the Neumann boundary conditions are listed
column2={'Boundary ID','Normal','Tangential'};
uitable('units','Normalized','Position',[0.845,0.10,0.13,0.44],...
    'Data',handles.dataNeu,...
    'ColumnName',column2,'RowName',[]);


%% ASSIGN DISPLACEMENT FUNCTION
% * Executes on button press in |assignDisp|;
% * Fills in the enforced load table for the boundaries selected in
% |listboxDir| with the string defined in |editDisp|, in the direction
% defined in the |popupmenu_NTDir|.

function assignDisp_Callback(hObject, eventdata, handles)
% hObject    handle to assignDisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
% Gets the selected boundaries in |listboxDir| and the direction of the
% load in |popupmenu_NTDir|
itemlist = get(handles.listboxDir,'Value'); % list of selected items in the listbox
nitems = length(itemlist); % number of selected items in the listbox
DispDirection = get(handles.popupmenu_NTDir,'Value');

%%
% Writes the displacemment definition in the column of |dataDir| associated
% to the direction the displacement is applied into.
for ii = 1:nitems
    crtitem = itemlist(ii);
    handles.dataDir{crtitem,DispDirection+1}=get(handles.editDisp,'string');
end

%%
% Updates the handles structure
guidata(hObject,handles);

%%
% Redraws the table where the Dirichlet boundary conditions are listed
column1={'Boundary ID','Normal','Tangential'};
uitable('units','Normalized','Position',[0.475,0.10,0.13,0.44],...
    'Data',handles.dataDir,...
    'ColumnName',column1,'RowName',[]);


%% RESET NEUMANN/DIRICHLET BOUNDARY CONDITION FUNCTIONS

%%
% *Reset Neumann*
%
% * Executes on button press in |resetNeumann|;
% * Substitutes all previous definitions in |dataNeu| by |NaN|;
% * Redraws the force boundary condition table.

function resetNeumann_Callback(hObject, eventdata, handles)
% hObject    handle to resetNeumann (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
% Substitutes all previous definitions in |dataNeu| by |NaN|
for i=1:length(handles.edgesNeumann)
    handles.dataNeu{i,2}='NaN';
    handles.dataNeu{i,3}='NaN';
end

%%
% Updates the handles structure
guidata(hObject,handles);

%%
% Redraws the table where the Neumann boundary conditions are listed
column2={'Boundary ID','Normal','Tangential'};
uitable('units','Normalized','Position',[0.845,0.10,0.13,0.44],...
    'Data',handles.dataNeu,...
    'ColumnName',column2,'RowName',[]);

%%
% *Reset Dirichlet*
%
% * Executes on button press in |resetDirichlet|;
% * Substitutes all previous definitions in |dataDir| by |NaN|;
% * Redraws the displacement boundary condition table.

function resetDirichlet_Callback(hObject, eventdata, handles)
% hObject    handle to resetDirichlet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
% Substitutes all previous definitions in |dataDir| by |NaN|
for i=1:length(handles.edgesDirichlet)
    handles.dataDir{i,2}='NaN';
    handles.dataDir{i,3}='NaN';
end

%%
% Updates the handles structure
guidata(hObject,handles);

%%
% Redraws the table where the Dirichlet boundary conditions are listed
column1={'Boundary ID','Normal','Tangential'};
uitable('units','Normalized','Position',[0.475,0.10,0.13,0.44],...
    'Data',handles.dataDir,...
    'ColumnName',column1,'RowName',[]);


%% PREVIOUS FUNCTION
% * Executes on button press in |previous|;
% * Just closes the current GUI and launches the previous one. All changes
% made in the current GUI are lost.

function previous_Callback(hObject, eventdata, handles)
% hObject    handle to previous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(StructTriBC2);

% Trying to close Visualize if it was opened
try
    close('Visualize');
catch
end

StructTriBC1;


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

VisualizeTri;


%% GUI ENTITIES GENERATION CODE
% * Automatically generated code for the buttons and menus;
% * Some (not-so-sound) checks are performed on the data inserted by the
% user in |editDisp| and |editNeu| just to make sure they are
% mathematically legible. 

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
    set(hObject,'string','');
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



