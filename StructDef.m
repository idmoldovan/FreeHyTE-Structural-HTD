%% PREDEFINED FUNCTIONS
% No changes were made to the following functions, predefined by GUIDE

function varargout = StructDef(varargin)
% STRUCTDEF MATLAB code for StructDef.fig
%      STRUCTDEF, by itself, creates a new STRUCTDEF or raises the existing
%      singleton*.
%
%      H = STRUCTDEF returns the handle to a new STRUCTDEF or the handle to
%      the existing singleton*.
%
%      STRUCTDEF('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STRUCTDEF.M with the given input arguments.
%
%      STRUCTDEF('Property','Value',...) creates a new STRUCTDEF or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before StructDef_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to StructDef_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help StructDef

% Last Modified by GUIDE v2.5 11-May-2016 10:58:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @StructDef_OpeningFcn, ...
                   'gui_OutputFcn',  @StructDef_OutputFcn, ...
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
% UIWAIT makes StructDef wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = StructDef_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%% OPENING FUNCTION
% * Executes just before StructDef is made visible;
% * Reads the previous data in the local |mat| file and fills in the
% fields.

function StructDef_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to StructDef (see VARARGIN)

%%
% If no varargin was defined, it is set to zero and saved to handles
if isempty(varargin)
    handles.varargin = 0;
else
    handles.varargin = varargin{1};
end

%%
% % Splash screen sequence
% if handles.varargin == 0
%     timer = tic;
%     ShowSplashForSecs = 1;
%     
%     s = SplashScreen( 'Splashscreen', 'Splash.png', ...
%         'ProgressBar', 'on', ...
%         'ProgressPosition', 5, ...
%         'ProgressRatio', 0.4 );
%     s.addText( 520, 80, 'Structural HTD', 'FontSize', 35, 'Color', [0 0 0.6] );
%     s.addText( 520, 120, 'v1.0', 'FontSize', 25, 'Color', [0.2 0.2 0.5] );
%     s.addText( 365, 270, 'Loading...', 'FontSize', 20, 'Color', 'white' );
%     
%     ElapsedTime = toc(timer);
%     if ElapsedTime < ShowSplashForSecs
%         pause(ShowSplashForSecs - ElapsedTime);
%     end
%     delete(s);
% end

%%
% Setting warnings to off. These warnings are caused by missing files and
% variables before the local |mat| files are written for the first time and
% by the possibility that the problem is purely Neumann or Dirichlet. The
% warnings are re-activated after the successful execution, at the end of
% the |main.m| function.
warning('off','MATLAB:DELETE:FileNotFound');
warning('off','MATLAB:load:variableNotFound');

%%
% Creates the |handles| structure
% Choose default command line output for |StructDef|
handles.output = hObject;

%% 
% Getting the current working folder (in R2012 may be different from the
% folder you started the app from!)
handles.WorkingFolder = pwd;

%%
% If there exists a local |mat| file, it loads its key elements to fill in
% the fields with data taken from the previous run.
if exist('./StructDef.mat','file')
    load('StructDef','L','B','Nx','Ny','EdgesOrder','LoopsOrder',...
        'NumberGaussPoints','MeshOption','Young','Poisson','PlaneState');
 
    %%
    % In rare situations, no values are stored for L, B, Nx, Ny. If this is
    % the case, it reads NaN and thus ChangesQ always results 1. To avoid
    % this, when NaN is red, it is substituted by zero.
    if isnan(L)
        L=0; B=0; Nx=0; Ny=0;
    end
    
    %%
    % Filling in the fields with the data from the previous iteration
    set(handles.edit_DimX,'String',sprintf('%d',L));
    set(handles.edit_DimY,'String',sprintf('%d',B));
    set(handles.edit_NLoopX,'String',sprintf('%d',Nx));
    set(handles.edit_NLoopY,'String',sprintf('%d',Ny));
    set(handles.edit_OrderEdge,'String',sprintf('%d',EdgesOrder));
    set(handles.edit_OrderLoop,'String',sprintf('%d',LoopsOrder));
    set(handles.edit_NGP,'String',sprintf('%d',NumberGaussPoints));
    set(handles.edit_Young,'String',sprintf('%g',Young));
    set(handles.edit_Poisson,'String',sprintf('%g',Poisson));
    set(handles.popupmenu_mesh,'Value',MeshOption);
    set(handles.popupmenu_plane,'Value',PlaneState);

    %%
    % If |MeshOption = 2|, that is, the mesh generation is automatic, it
    % makes no sense to edit the fields associated to the regular mesh
    % generator. They become inactive.
    if MeshOption == 1
        set(handles.edit_DimX, 'enable', 'on');
        set(handles.edit_DimY, 'enable', 'on');
        set(handles.edit_NLoopX, 'enable', 'on');
        set(handles.edit_NLoopY, 'enable', 'on');
    else
        set(handles.edit_DimX, 'enable', 'off');
        set(handles.edit_DimY, 'enable', 'off');
        set(handles.edit_NLoopX, 'enable', 'off');
        set(handles.edit_NLoopY, 'enable', 'off');
    end 
end

%% 
% Returns control to the user.
% Update handles structure
guidata(hObject, handles);





%% NEXT BUTTON
% * Executes on button press in |pushbutton_next|;
% * It recovers all data provided by the user in the GUI fields;
% * If the mesh is regular, it computes the |edgesArray| list (consisting of
% the exterior edges of the srtucture);
% * If the mesh is formed by triangular elements, it launches the |pdetool|
% to generate the mesh;
% * In either case, after the mesh generation is complete, it checks for
% relevant (i.e. mesh) changes as compared to the previous |mat| file. If
% such changes exist, it deletes the |mat| file corresponding to the next
% GUIs to force their definition from scratch. If no mesh changes were
% detected, the next GUI will load its previous version;
% * If the user asked for the model to be saved, it saves the information
% regarding this GUI in the specified file and folder.
function pushbutton_next_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
% Reading the information provided by the user in the GUI fields.
EdgesOrder = str2double(get(handles.edit_OrderEdge,'String'));
LoopsOrder = str2double(get(handles.edit_OrderLoop,'String'));
NumberGaussPoints = str2double(get(handles.edit_NGP,'String'));
MeshOption = get(handles.popupmenu_mesh,'Value');
Young = str2double(get(handles.edit_Young,'String'));
Poisson = str2double(get(handles.edit_Poisson,'String'));
PlaneState = get(handles.popupmenu_plane,'Value');

%%
% If the user asked for the model to be saved, it stores the path and the
% file name.
if isfield(handles,'DirName')
    DirName = handles.DirName;
    FileName = handles.FileName;
else
    DirName = '';
    FileName = '';
end

%% 
% *Procedure for the regular rectangular mesh*
if MeshOption == 1 
    %%
    % Reading mesh information
    L = str2double(get(handles.edit_DimX,'String'));
    B = str2double(get(handles.edit_DimY,'String'));
    Nx = str2double(get(handles.edit_NLoopX,'String'));
    Ny = str2double(get(handles.edit_NLoopY,'String'));
    %%
    % Computing the |edgesArray| list
    L1 = linspace(1,Ny,Ny);
    L2 = linspace((Nx*Ny)+1,((Nx*Ny)+1)+(Ny-1),Ny);
    L3 = linspace(((Nx*Ny)+1)+Ny,((Nx*Ny)+1)+Ny+((Ny+1)*(Nx-1)),Nx);
    L4 = linspace(((Nx*Ny)+1)+2*Ny,((Nx*Ny)+1)+Ny+((Ny+1)*(Nx-1))+Ny,Nx);
    edgesArray = sort(cat(1,L1',L2',L3',L4'));
    %%
    % |p| and |t| variables are allocated dummy values to avoid errors
    % related to the comparison between |mat| files corresponding to
    % regular (rectangular) and triangular meshes.
    p = 0; t = 0; 
    %%
    % Check if critical changes were made in the current GUI session. Note
    % that critical changes only refer to the mesh, not necessarily to the
    % geometry of the structure.
    if exist('./StructDef.mat','file')
        save('ToDetectChanges','MeshOption','Nx','Ny','p','t'); % temporary file to detect critical changes
        % Fields whose change triggers the change of the mesh
        CriticalFields = {'MeshOption','Nx','Ny','p','t'};
        % Checks the old |StructDef.m| and the new |ToDetectChanges| file
        % for changes in the critical fields. ChangesQ = 1 if there are 
        % changes, 0 otherwise.
        ChangesQ = any(~isfield(comp_struct(load('StructDef.mat'),...
            load('ToDetectChanges.mat')),CriticalFields));
        % deleting the auxiliary file
        delete('ToDetectChanges.mat');
    else
        % If there is no |StructDef| file in the first place, it sets
        % ChangesQ to 1.
        ChangesQ = 1;
    end
    
    %%        
    % Saving the workspace to the local |mat| file. If requested by the 
    % user, the GUI information is also saved to the file specified by the
    % user.
    save('StructDef','L','B','Nx','Ny','EdgesOrder','LoopsOrder',...
        'NumberGaussPoints','MeshOption','Young','Poisson','PlaneState',...
        'p','t','edgesArray','DirName','FileName','ChangesQ');
    
    % If the folder where it looks for pre-load files is different from the
    % current folder, it creates a copy of StructDef in the former
    if  ~strcmp(pwd,handles.WorkingFolder)
        save(fullfile(handles.WorkingFolder,'StructDef'),'L','B','Nx','Ny',...
            'EdgesOrder','LoopsOrder',...
            'NumberGaussPoints','MeshOption','Young','Poisson','PlaneState',...
            'p','t','edgesArray','DirName','FileName','ChangesQ');
    end
    
    if ~isempty(DirName)
       save(fullfile(DirName,FileName),...
            'L','B','Nx','Ny','EdgesOrder','LoopsOrder',...
            'NumberGaussPoints','MeshOption','Young','Poisson','PlaneState',...
            'p','t','edgesArray','ChangesQ'); 
    end
        
    %%
    % Preparing to exit.
    close(handles.figure1);
        
    %% 
    % If there are relevant changes, it physically deletes the |mat| files
    % associated to the following GUIs.
    if ChangesQ
        delete('StructRegBC1.mat','StructRegBC2.mat');
    end
    StructRegBC1; 

%% 
% *Procedure for the automatic triangular mesh*    
else 
    %%
    % Reading mesh information. This data is useless for the automatic mesh
    % generation. It is only stored to fill in the corresponding fields in
    % a future run.
    L = str2double(get(handles.edit_DimX,'String'));
    B = str2double(get(handles.edit_DimY,'String'));
    Nx = str2double(get(handles.edit_NLoopX,'String'));
    Ny = str2double(get(handles.edit_NLoopY,'String'));
    %%
    % |edgesArray| variable is allocated a dummy value to avoid errors
    % related to the comparison between |mat| files corresponding to
    % regular (rectangular) and triangular meshes.
    edgesArray = 0; 
    %%
    % Launching |pdetool| to define the mesh
    pdewindow = pdeinit; 
    %% 
    % Stopping the execution until the |pdetool| window is closed 
    waitfor(pdewindow);
    %%
    % This is a basic check to confirm that the user saved the mesh 
    % information. A new mesh need not be created if the corresponding
    % information already exist in |StructDef.mat|, so the program checks
    % if a new mesh was defined or if |StructDef.mat| exists. Of course,
    % the check fails if |StructDef.mat| exists, but does not contain
    % relevant mesh information (the program exists with an error).
    BaseVars = evalin('base','whos'); % collects variable info from base
    % if the mesh variables are not found in base and an old definition
    % does not exist in StructDef.mat
    if (~ismember('p',[BaseVars(:).name]) || ~ismember('t',[BaseVars(:).name])) && ...
            (~exist('./StructDef.mat','file'))
        errordlg('No mesh information was exported or variable names were changed. Please press "Next" again and export the mesh info as p, e and t.','Invalid input','modal');
        uicontrol(hObject);
        return;
    end
    
    %%
    % if the mesh variables are found in base, it loads them...
    if ismember('p',[BaseVars(:).name]) && ismember('t',[BaseVars(:).name])
        p = evalin('base','p');
        t = evalin('base','t');
    %%
    % ... otherwise, it loads them from the StructDef.mat
    else
        load('StructDef','p','t');
    end
    
    %%
    % Check if critical changes were made in the current GUI session. Note
    % that critical changes only refer to the mesh, not necessarily to the
    % geometry of the structure.
    if exist('./StructDef.mat','file')
        save('ToDetectChanges','MeshOption','Nx','Ny','p','t'); % temporary file to detect critical changes
        % Fields whose change triggers the change of the mesh
        CriticalFields = {'MeshOption','Nx','Ny','p','t'};
        % Checks the old StructDef and the new ToDetectChanges for changes in
        % the critical fields. ChangesQ = 1 if there are changes, 0 otherwise
        ChangesQ = any(~isfield(comp_struct(load('StructDef.mat'),...
            load('ToDetectChanges.mat')),CriticalFields));
        % deleting the auxiliary file
        delete('ToDetectChanges.mat');
    else
        % If there is no |StructDef| file in the first place, it sets
        % ChangesQ to 1.
        ChangesQ = 1;
    end
    
        %%        
    % Saving the workspace to the local |mat| file. If requested by the 
    % user, the GUI information is also saved to the file specified by the
    % user.
    save('StructDef','L','B','Nx','Ny','EdgesOrder','LoopsOrder',...
        'NumberGaussPoints','MeshOption','Young','Poisson','PlaneState',...
        'p','t','edgesArray','DirName','FileName','ChangesQ');
    
    % If the folder where it looks for pre-load files is different from the
    % current folder, it creates a copy of StructDef in the former
    if  ~strcmp(pwd,handles.WorkingFolder)
        save(fullfile(handles.WorkingFolder,'StructDef'),'L','B','Nx','Ny',...
            'EdgesOrder','LoopsOrder',...
            'NumberGaussPoints','MeshOption','Young','Poisson','PlaneState',...
            'p','t','edgesArray','DirName','FileName','ChangesQ');
    end
    
    % if requested by the user, saves to file
    if ~isempty(DirName)
        save(fullfile(DirName,FileName),...
            'L','B','Nx','Ny','EdgesOrder','LoopsOrder',...
            'NumberGaussPoints','MeshOption','Young','Poisson','PlaneState',...
            'p','t','edgesArray','ChangesQ');
    end
    
    %%
    % Preparing to exit.
    close(handles.figure1);
    
    %% 
    % If there are relevant changes, it physically deletes the |mat| files
    % associated to the following GUIs.
    if ChangesQ
        delete('StructTriBC1.mat','StructTriBC2.mat');
    end
    StructTriBC1; 
end


%% RESET BUTTON
% * Executes on button press in |pushbutton_reset|;
% * It deletes all entries of the edit fields and resets the pop-up menus;
% * It also deletes the save information so that the save button action is
% reversed.
function pushbutton_reset_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
% Removing all field information
set(handles.edit_DimX,'String','');
set(handles.edit_DimY,'String','');
set(handles.edit_NLoopX,'String','');
set(handles.edit_NLoopY,'String','');
set(handles.edit_OrderEdge,'String','');
set(handles.edit_OrderLoop,'String','');
set(handles.edit_NGP,'String','');
set(handles.edit_Young,'String','');
set(handles.edit_Poisson,'String','');

%%
% Reseting |edit_Path|, the pop-up menus, and the properties of the mesh
% definition fields
set(handles.edit_Path,'String',...
    'THE MODEL WILL NOT BE SAVED!  Please press Save if you wish to save the model.',...
    'BackgroundColor',[1.0 0.694 0.392]);
set(handles.popupmenu_plane,'Value',1);
set(handles.popupmenu_mesh,'Value',1);    

set(handles.edit_DimX, 'enable', 'on');
set(handles.edit_DimY, 'enable', 'on');
set(handles.edit_NLoopX, 'enable', 'on');
set(handles.edit_NLoopY, 'enable', 'on');

handles.FileName = '';
handles.DirName = '';

%%
% Updating the |handles| structure
guidata(hObject, handles);


%% SAVE BUTTON
% * Executes on button press in |Save|;
% * Reads the path and filename provided by the user to save the model;
% * Generates the |FileName| and |DirName| members in the |handles|
% structure to store this information and to make it available to other
% functions.
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% 
% Reads the |handles| structure
handles=guidata(hObject);

%%
% Gets the path and filename to save the model
[FileName,DirName] = uiputfile('*.mat','Save as');

if FileName
    %%
    % Registers the path to the save file in the |edit_Path| box
    set(handles.edit_Path,'string',fullfile(DirName,FileName),'BackgroundColor',[0 1 0]);
    
    %%
    % Generates the |FileName| and |DirName| members in the |handles|
    % structure to store this information and to make it available to other
    % functions.
    handles.FileName = FileName;
    handles.DirName = DirName;
    
    %%
    % If user cancelled the saving process and no save file was given
else
    handles.FileName = '';
    handles.DirName = '';
    
end

%%
% Updating the |handles| structure
guidata(hObject, handles);


%% LOAD BUTTON
% * Executes on button press in |Load|;
% * Loads all relevant variables in the file appointed by the user into the 
% workspace. If the definition of the model was not completed in the
% previous session, some of the variables cannot be loaded (the
% corresponding warning is suppressed);
% * Stores the required variables into the StructDef |mat| file. This data
% must always be available;
% * Verifies if the data generated by the second and third GUIs are
% available. If they are, stores them into the local |mat| files, otherwise
% _deletes_ the |mat| files, to force their reinitialization;
% * Deletes the |FileName| and |DirName| variables to avoid overwriting of the
% loaded file;
% * Refreshes the interface, filling in the fields with the newly loaded
% values.
function Load_Callback(hObject, eventdata, handles)
% hObject    handle to Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% 
% Reads the |handles| structure
handles=guidata(hObject);

%%
% Loads the |mat| file indicated by the user
[FileName,DirName] = uigetfile('*.mat','File to load');

%%
% Deletes the |FileName| and |DirName| variables to avoid overwriting;
handles.FileName = '';
handles.DirName = '';

%%
% Loading all relevant variables into the current workspace
load(fullfile(DirName,FileName),'L','B','Nx','Ny','EdgesOrder','LoopsOrder',...
    'NumberGaussPoints','MeshOption','Young','Poisson','PlaneState','p','t',...
    'ChangesQ','edgesArray','data','dataDir','dataNeu','edgesDirichlet',...
    'edgesNeumann');

%%
% Saving the local |StructDef.mat| file 
save('StructDef','L','B','Nx','Ny','EdgesOrder','LoopsOrder',...
    'NumberGaussPoints','MeshOption','Young','Poisson','PlaneState','p','t',...
    'edgesArray','ChangesQ');

% If the folder where it looks for pre-load files is different from the
% current folder, it creates a copy of StructDef in the former
if  ~strcmp(pwd,handles.WorkingFolder)
    save(fullfile(handles.WorkingFolder,'StructDef'),'L','B','Nx','Ny',...
        'EdgesOrder','LoopsOrder',...
        'NumberGaussPoints','MeshOption','Young','Poisson','PlaneState',...
        'p','t','edgesArray','DirName','FileName','ChangesQ');
end

%%
% If 'p' exists (is not zero), updates 'p' and 't' in the base workspace
if any(any(p~=0))
    assignin('base','p',p);
    assignin('base','t',t);
end

%%
% Depending on the kind of mesh that is selected ...
if MeshOption == 1 % regular rectangular mesh
    %%
    % ... checks if the variable created by the second GUI exists and if it
    % doesn't, it deletes the local |mat| file to force reinitialization...
    if ~exist('data','var')
        delete('StructRegBC1.mat');
        %%
        % ... or it stores (and overwrites) the |mat| file if the variable
        % exists.
    else
        save('StructRegBC1','data');
    end
    
    %%
    % ... checks if a variable created by the third GUI exists and if it
    % doesn't, it deletes the local |mat| file to force reinitialization...
    if ~exist('dataDir','var') && ~exist('dataNeu','var')
        delete('StructRegBC2.mat');
        %%
        % ... or it stores (and overwrites) the |mat| file if the variable
        % exists.
    else
        delete('StructRegBC2.mat');
        dummy = 0;
        save('StructRegBC2','dummy');
        if exist('dataDir','var')
            save('StructRegBC2','-append','edgesDirichlet','dataDir');
        end
        if exist('dataNeu','var')
            save('StructRegBC2','-append','edgesNeumann','dataNeu');
        end       
    end
    
%%
% Same operation for the irregular mesh file.
else
    if ~exist('data','var')
        delete('StructTriBC1.mat');
    else
        save('StructTriBC1','data');
    end
    
    if ~exist('dataDir','var') && ~exist('dataNeu','var')
        delete('StructTriBC2.mat');
    else
        delete('StructTriBC2.mat');
        dummy = 0;
        save('StructTriBC2','dummy');
        if exist('dataDir','var')
            save('StructTriBC2','-append','edgesDirichlet','dataDir');
        end
        if exist('dataNeu','var')
            save('StructTriBC2','-append','edgesNeumann','dataNeu');
        end
    end
end
 
%%
% Refreshes the interface, writing the loaded values into its fields...
set(handles.edit_DimX,'String',sprintf('%d',L));
set(handles.edit_DimY,'String',sprintf('%d',B));
set(handles.edit_NLoopX,'String',sprintf('%d',Nx));
set(handles.edit_NLoopY,'String',sprintf('%d',Ny));
set(handles.edit_OrderEdge,'String',sprintf('%d',EdgesOrder));
set(handles.edit_OrderLoop,'String',sprintf('%d',LoopsOrder));
set(handles.edit_NGP,'String',sprintf('%d',NumberGaussPoints));
set(handles.edit_Young,'String',sprintf('%g',Young));
set(handles.edit_Poisson,'String',sprintf('%g',Poisson));
set(handles.popupmenu_mesh,'Value',MeshOption);
set(handles.popupmenu_plane,'Value',PlaneState);
set(handles.edit_Path,'string',...
    'THE MODEL WILL NOT BE SAVED!  Please press Save if you wish to save the model.',...
    'BackgroundColor',[1.0 0.694 0.392]);

%%
% ... and activates or inactivates the regular mesh fields according to the
% type of mesh in the loaded model.
if MeshOption == 1
    set(handles.edit_DimX, 'enable', 'on');
    set(handles.edit_DimY, 'enable', 'on');
    set(handles.edit_NLoopX, 'enable', 'on');
    set(handles.edit_NLoopY, 'enable', 'on');
else
    set(handles.edit_DimX, 'enable', 'off');
    set(handles.edit_DimY, 'enable', 'off');
    set(handles.edit_NLoopX, 'enable', 'off');
    set(handles.edit_NLoopY, 'enable', 'off');
end

%%
% Updates the |handles| structure
guidata(hObject, handles);


%% CLEAR BUTTON
% * Executes on button press in |Clear|;
% * Deletes the path and name of the save file and reinitializes the
% |edit_Path| field.
function Clear_Callback(hObject, eventdata, handles)
% hObject    handle to Clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% 
% Reads the |handles| structure
handles=guidata(hObject);

%%
% Deletes the |FileName| and |DirName| variables 
handles.FileName = '';
handles.DirName = '';
set(handles.edit_Path,'string',...
    'THE MODEL WILL NOT BE SAVED!  Please press Save if you wish to save the model.',...
    'BackgroundColor',[1.0 0.694 0.392]);

%%
% Updates the |handles| structure
guidata(hObject, handles);


%% POP-UP MENUS
% *Mesh pop-up menu*
% * Executes on selection change in |popupmenu_mesh|.
% * Activates/deactivates the regular mesh definition fields according to
% the type of mesh that the user chooses.
function popupmenu_mesh_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_mesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_mesh contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_mesh
if get(handles.popupmenu_mesh,'Value') == 1
    set(handles.edit_DimX, 'enable', 'on');
    set(handles.edit_DimY, 'enable', 'on');
    set(handles.edit_NLoopX, 'enable', 'on');
    set(handles.edit_NLoopY, 'enable', 'on');
else
   set(handles.edit_DimX, 'enable', 'off');
   set(handles.edit_DimY, 'enable', 'off');
   set(handles.edit_NLoopX, 'enable', 'off');
   set(handles.edit_NLoopY, 'enable', 'off');
end


% --- Executes during object creation, after setting all properties.
function popupmenu_mesh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_mesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% 
% *Plane pop-up menu*
% * Executes on selection change in popupmenu_plane;
% * Values 1 and 2 correspond to plane stress and plane strain,
% respectively.
function popupmenu_plane_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_plane contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_plane

% --- Executes during object creation, after setting all properties.
function popupmenu_plane_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end







%% EDIT FIELDS
% Field box editing functions. Most Callbacks check for the validity of the
% data.

function edit_DimX_Callback(hObject, eventdata, handles)
% hObject    handle to edit_DimX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_DimX as text
%        str2double(get(hObject,'String')) returns contents of edit_DimX as a double
L = str2double(get(hObject,'String'));
if isnan(L) || ~isreal(L) || isequal(L,0) || L < 0
    set(hObject,'String','');
    errordlg('You must enter a real positive value','Invalid input','modal');
    uicontrol(hObject);
    return;
end
% --- Executes during object creation, after setting all properties.
function edit_DimX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_DimX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_DimY_Callback(hObject, eventdata, handles)
% hObject    handle to edit_DimY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_DimY as text
%        str2double(get(hObject,'String')) returns contents of edit_DimY as a double
B = str2double(get(hObject,'String'));
if isnan(B) || ~isreal(B) || isequal(B,0) || B < 0
    set(hObject,'String','');
    errordlg('You must enter a real positive value','Invalid input','modal');
    uicontrol(hObject);
    return;
end
% --- Executes during object creation, after setting all properties.
function edit_DimY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_DimY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_NLoopX_Callback(hObject, eventdata, handles)
% hObject    handle to edit_NLoopX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_NLoopX as text
%        str2double(get(hObject,'String')) returns contents of edit_NLoopX as a double
Nx = str2double(get(hObject,'String'));
if isnan(Nx) || ~isreal(Nx) || logical(abs(round(Nx)-Nx)<eps)==0 || isequal(Nx,0) || Nx < 0
    set(hObject,'String','');
    errordlg('You must enter a positive integer value','Invalid input','modal');
    uicontrol(hObject);
    return;
end
% --- Executes during object creation, after setting all properties.
function edit_NLoopX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_NLoopX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_NLoopY_Callback(hObject, eventdata, handles)
% hObject    handle to edit_NLoopY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_NLoopY as text
%        str2double(get(hObject,'String')) returns contents of edit_NLoopY as a double
Ny = str2double(get(hObject,'String'));
if isnan(Ny) || ~isreal(Ny) || logical(abs(round(Ny)-Ny)<eps)==0 || isequal(Ny,0) || Ny < 0
    set(hObject,'String','');
    errordlg('You must enter a positive integer value','Invalid input','modal');
    uicontrol(hObject);
    return;
end
% --- Executes during object creation, after setting all properties.
function edit_NLoopY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_NLoopY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_OrderEdge_Callback(hObject, eventdata, handles)
% hObject    handle to edit_OrderEdge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_OrderEdge as text
%        str2double(get(hObject,'String')) returns contents of edit_OrderEdge as a double
EdgesOrder = str2double(get(handles.edit_OrderEdge,'String'));
if isnan(EdgesOrder) || ~isreal(EdgesOrder) ||... 
        logical(abs(round(EdgesOrder)-EdgesOrder)<eps)==0 || EdgesOrder < 0
    set(hObject,'String','');
    errordlg('You must enter either 0 or a natural number','Invalid input','modal');
    uicontrol(hObject);
    return;
end
% --- Executes during object creation, after setting all properties.
function edit_OrderEdge_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_OrderEdge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_OrderLoop_Callback(hObject, eventdata, handles)
% hObject    handle to edit_OrderLoop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_OrderLoop as text
%        str2double(get(hObject,'String')) returns contents of edit_OrderLoop as a double
LoopsOrderC = str2double(get(handles.edit_OrderLoop,'String'));
if isnan(LoopsOrderC) || ~isreal(LoopsOrderC) || ...
        logical(abs(round(LoopsOrderC)-LoopsOrderC)<eps)==0 || LoopsOrderC < 0
    set(hObject,'String','');
    errordlg('You must enter either 0 or a natural number','Invalid input','modal');
    uicontrol(hObject);
    return;
end
% --- Executes during object creation, after setting all properties.
function edit_OrderLoop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_OrderLoop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_NGP_Callback(hObject, eventdata, handles)
% hObject    handle to edit_NGP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_NGP as text
%        str2double(get(hObject,'String')) returns contents of edit_NGP as a double
NumberGaussPoints = str2double(get(handles.edit_NGP,'String'));
if isnan(NumberGaussPoints) || ~isreal(NumberGaussPoints) || ...
        logical(abs(round(NumberGaussPoints)-NumberGaussPoints)<eps)==0 ||... 
        isequal(NumberGaussPoints,0) || NumberGaussPoints < 0
    set(hObject,'String','');
    errordlg('You must enter a positive integer value','Invalid input','modal');
    uicontrol(hObject);
    return;
end
% --- Executes during object creation, after setting all properties.
function edit_NGP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_NGP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Poisson_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Poisson (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Poisson as text
%        str2double(get(hObject,'String')) returns contents of edit_Poisson as a double
Poisson = str2double(get(handles.edit_Poisson,'String'));
if isnan(Poisson) || ~isreal(Poisson) || Poisson <= 0 || Poisson >= 0.5
    set(hObject,'String','');
    errordlg('You must enter a positive real number smaller than 0.5','Invalid input','modal');
    uicontrol(hObject);
    return;
end
% --- Executes during object creation, after setting all properties.
function edit_Poisson_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Poisson (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Young_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Young (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Young as text
%        str2double(get(hObject,'String')) returns contents of edit_Young as a double
Young = str2double(get(handles.edit_Young,'String'));
if isnan(Young) || ~isreal(Young) || Young <= 0
    set(hObject,'String','');
    errordlg('You must enter a positive real number','Invalid input','modal');
    uicontrol(hObject);
    return;
end
% --- Executes during object creation, after setting all properties.
function edit_Young_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Young (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_Path_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_Path as text
%        str2double(get(hObject,'String')) returns contents of edit_Path as a double

% --- Executes during object creation, after setting all properties.
function edit_Path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




%% TEXT FIELDS
% Text fields' editing functions. The respective |ButtonDownFcn| defines
% the help of that field, accessed through right clicking.

% --- Executes during object creation, after setting all properties.
function text_DimX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_DimX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function text_DimY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_DimY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function text_NLoopX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_NLoopX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function text_NLoopY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_NLoopY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function text_OrderEdge_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_OrderEdge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function text_OrderLoop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_OrderLoop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function text_NGP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_NGP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function text_SaveLoad_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_SaveLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
