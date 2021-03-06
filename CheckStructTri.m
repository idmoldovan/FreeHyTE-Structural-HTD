function varargout = CheckStructTri(varargin)
% CHECKSTRUCTTRI MATLAB code for CheckStructTri.fig
%      CHECKSTRUCTTRI, by itself, creates a new CHECKSTRUCTTRI or raises the existing
%      singleton*.
%
%      H = CHECKSTRUCTTRI returns the handle to a new CHECKSTRUCTTRI or the handle to
%      the existing singleton*.
%
%      CHECKSTRUCTTRI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHECKSTRUCTTRI.M with the given input arguments.
%
%      CHECKSTRUCTTRI('Property','Value',...) creates a new CHECKSTRUCTTRI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CheckStructTri_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CheckStructTri_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CheckStructTri

% Last Modified by GUIDE v2.5 13-Jul-2016 22:00:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CheckStructTri_OpeningFcn, ...
                   'gui_OutputFcn',  @CheckStructTri_OutputFcn, ...
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


% --- Executes just before CheckStructTri is made visible.
function CheckStructTri_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CheckStructTri (see VARARGIN)

% Choose default command line output for CheckStructTri
handles.output = hObject;

% Launch input processor
[~,Nodes,Edges,Loops,~,~] = InputProcTri;

% prepare to plot
xmesh = [Nodes(Edges.nini(:),1) Nodes(Edges.nfin(:),1)];
ymesh = [Nodes(Edges.nini(:),2) Nodes(Edges.nfin(:),2)];
% get delta to scale the insertion points of the text
delta = sqrt(max(max(xmesh))^2+max(max(ymesh))^2); 

for ii=1:length(Edges.type)
    line(xmesh(ii,:),ymesh(ii,:),'LineWidth',1,'color','b');
end
daspect([1 1 1]);
axis off ; 

% Plotting the edges
for i = 1:length(Edges.type)
    EX = [Nodes(Edges.nini(i),1); Nodes(Edges.nfin(i),1)];
    EY = [Nodes(Edges.nini(i),2); Nodes(Edges.nfin(i),2)];
    pos = [sum(EX)/2+0.0*delta,sum(EY)/2+0.0*delta] ;
    text(pos(1),pos(2),int2str(i),'fontsize',8, ...
        'BackgroundColor','w','fontweight','bold','color','r');
end

% Plotting the elements
for i = 1:length(Loops.area)
    C = polygonCentroid(Nodes(Loops.nodes(i,:),:));
    pos = C ;
    text(pos(1),pos(2),int2str(i),'fontsize',8, 'fontweight','bold',...
        'BackgroundColor','w','color','b');
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CheckStructTri wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CheckStructTri_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(CheckStructTri); % closing the GUI window
pause(0.05); % to actually close it
MainTri;

% --- Executes on button press in previous.
function previous_Callback(hObject, eventdata, handles)
% hObject    handle to previous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(CheckStructTri);
StructTriBC2;

% --- Executes on button press in update.
function update_Callback(hObject, eventdata, handles)
% hObject    handle to update (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% clearing the display
cla(gcf) ;

% getting the ViewOption
ViewOption = get(handles.view,'Value');

% Launch input processor
[~,Nodes,Edges,Loops,BConds,~] = InputProcTri;

% prepare to plot
xmesh = [Nodes(Edges.nini(:),1) Nodes(Edges.nfin(:),1)];
ymesh = [Nodes(Edges.nini(:),2) Nodes(Edges.nfin(:),2)];
% get delta to scale the insertion points of the text
delta = sqrt(max(max(xmesh))^2+max(max(ymesh))^2);

if ViewOption==1 % just the mesh
    
    for ii=1:length(Edges.type)
        line(xmesh(ii,:),ymesh(ii,:),'LineWidth',1,'color','b');
    end
    daspect([1 1 1]);
    
    % Plotting the edges
    for i = 1:length(Edges.type)
        EX = [Nodes(Edges.nini(i),1); Nodes(Edges.nfin(i),1)];
        EY = [Nodes(Edges.nini(i),2); Nodes(Edges.nfin(i),2)];
        pos = [sum(EX)/2+0.0*delta,sum(EY)/2+0.0*delta] ;
        text(pos(1),pos(2),int2str(i),'fontsize',8, ...
            'BackgroundColor','w','fontweight','bold','color','r');
    end
    
    % Plotting the elements
    for i = 1:length(Loops.area)
        C = polygonCentroid(Nodes(Loops.nodes(i,:),:));
        pos = C ;
        text(pos(1),pos(2),int2str(i),'fontsize',8, 'fontweight','bold',...
            'BackgroundColor','w','color','b');
    end
    
elseif ViewOption==2 % just the BC, normal
    
    for ii=1:length(Edges.type)
        if strcmp(Edges.type(ii),'D') && ~Edges.lright(ii)
            line(xmesh(ii,:),ymesh(ii,:),'LineWidth',3,'color','k');
            diffx = xmesh(ii,2)-xmesh(ii,1);
            diffy = ymesh(ii,2)-ymesh(ii,1);
            posini = Nodes(Edges.nini(ii),:)+ [0.08*diffx 0.08*diffy];
            posend = Nodes(Edges.nfin(ii),:)- [0.08*diffx 0.08*diffy];
            text(posini(1),posini(2),num2str(BConds.Dirichlet{ii,1}(1)),...
                'HorizontalAlignment','center',...
                'BackgroundColor','w','fontsize',7,'color','k');
            text(posend(1),posend(2),num2str(BConds.Dirichlet{ii,1}(end)),...
                'HorizontalAlignment','center',...
                'BackgroundColor','w','fontsize',7,'color','k');
        elseif strcmp(Edges.type(ii),'D') && Edges.lright(ii)
            line(xmesh(ii,:),ymesh(ii,:),'LineWidth',1,'color','b');
        else
            line(xmesh(ii,:),ymesh(ii,:),'LineWidth',3,'color','r');
            diffx = xmesh(ii,2)-xmesh(ii,1);
            diffy = ymesh(ii,2)-ymesh(ii,1);
            posini = Nodes(Edges.nini(ii),:)+ [0.08*diffx 0.08*diffy];
            posend = Nodes(Edges.nfin(ii),:)- [0.08*diffx 0.08*diffy];
            text(posini(1),posini(2),num2str(BConds.Neumann{ii,1}(1)),...
                'HorizontalAlignment','center',...
                'BackgroundColor','w','fontsize',7,'color','r');
            text(posend(1),posend(2),num2str(BConds.Neumann{ii,1}(end)),...
                'HorizontalAlignment','center',...
                'BackgroundColor','w','fontsize',7,'color','r');
        end
    end
else % just the BC, tangential
    
    for ii=1:length(Edges.type)
        if strcmp(Edges.type(ii),'D') && ~Edges.lright(ii)
            line(xmesh(ii,:),ymesh(ii,:),'LineWidth',3,'color','k');
            diffx = xmesh(ii,2)-xmesh(ii,1);
            diffy = ymesh(ii,2)-ymesh(ii,1);
            posini = Nodes(Edges.nini(ii),:)+ [0.08*diffx 0.08*diffy];
            posend = Nodes(Edges.nfin(ii),:)- [0.08*diffx 0.08*diffy];
            text(posini(1),posini(2),num2str(BConds.Dirichlet{ii,2}(1)),...
                'HorizontalAlignment','center',...
                'BackgroundColor','w','fontsize',7,'color','k');
            text(posend(1),posend(2),num2str(BConds.Dirichlet{ii,2}(end)),...
                'HorizontalAlignment','center',...
                'BackgroundColor','w','fontsize',7,'color','k');
        elseif strcmp(Edges.type(ii),'D') && Edges.lright(ii)
            line(xmesh(ii,:),ymesh(ii,:),'LineWidth',1,'color','b');
        else
            line(xmesh(ii,:),ymesh(ii,:),'LineWidth',3,'color','r');
            diffx = xmesh(ii,2)-xmesh(ii,1);
            diffy = ymesh(ii,2)-ymesh(ii,1);
            posini = Nodes(Edges.nini(ii),:)+ [0.08*diffx 0.08*diffy];
            posend = Nodes(Edges.nfin(ii),:)- [0.08*diffx 0.08*diffy];
            text(posini(1),posini(2),num2str(BConds.Neumann{ii,2}(1)),...
                'HorizontalAlignment','center',...
                'BackgroundColor','w','fontsize',7,'color','r');
            text(posend(1),posend(2),num2str(BConds.Neumann{ii,2}(end)),...
                'HorizontalAlignment','center',...
                'BackgroundColor','w','fontsize',7,'color','r');
        end
    end
end



% --- Executes on selection change in view.
function view_Callback(hObject, eventdata, handles)
% hObject    handle to view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns view contents as cell array
%        contents{get(hObject,'Value')} returns selected item from view


% --- Executes during object creation, after setting all properties.
function view_CreateFcn(hObject, eventdata, handles)
% hObject    handle to view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
