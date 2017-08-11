function varargout = VisualizeReg(varargin)
% VISUALIZEREG MATLAB code for VisualizeReg.fig
%      VISUALIZEREG, by itself, creates a new VISUALIZEREG or raises the existing
%      singleton*.
%
%      H = VISUALIZEREG returns the handle to a new VISUALIZEREG or the handle to
%      the existing singleton*.
%
%      VISUALIZEREG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VISUALIZEREG.M with the given input arguments.
%
%      VISUALIZEREG('Property','Value',...) creates a new VISUALIZEREG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before VisualizeReg_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to VisualizeReg_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help VisualizeReg

% Last Modified by GUIDE v2.5 15-Jul-2016 10:18:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @VisualizeReg_OpeningFcn, ...
                   'gui_OutputFcn',  @VisualizeReg_OutputFcn, ...
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


% --- Executes just before VisualizeReg is made visible.
function VisualizeReg_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to VisualizeReg (see VARARGIN)

% Choose default command line output for VisualizeReg
handles.output = hObject;

% Create data structures
load('Visualize','nodes','edges_nodes','edges_loops','loops_nodes',...
     'loops_edges');

% Definitions of the data structures
Nodes = nodes;
Edges=struct('nini',edges_nodes(:,1),'nfin',edges_nodes(:,2),...
    'parametric',createLine(Nodes(edges_nodes(:,1),:),...
    Nodes(edges_nodes(:,2),:)),...
    'lleft',edges_loops(:,1),'lright',edges_loops(:,2),...
    'type',char(zeros(length(edges_nodes(:,1)),1)),'order',...
    zeros(length(edges_nodes(:,1)),1)); 
Loops=struct('nodes',loops_nodes,'edges',loops_edges,...
    'center',zeros(length(loops_nodes(:,1)),2),...
    'area',zeros(length(loops_nodes(:,1)),1));


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

% UIWAIT makes VisualizeReg wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = VisualizeReg_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
