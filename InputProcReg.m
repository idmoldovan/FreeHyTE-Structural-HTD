function [NGP, Nodes, Edges, Loops, BConds, NoDiv] = ...
    InputProcReg
% Input processing function for regular meshes

%% ------------------ User Input Zone ------------------

% Loads mesh information from previous interfaces 
load('StructRegBC1','nodes','edges_nodes','edges_loops','loops_nodes',...
     'loops_edges');

% Definitions of the data structures
Nodes = nodes;
Edges=struct('nini',edges_nodes(:,1),'nfin',edges_nodes(:,2),...
    'parametric',createLine(Nodes(edges_nodes(:,1),:),...
    Nodes(edges_nodes(:,2),:)),...
    'lleft',edges_loops(:,1),'lright',edges_loops(:,2),...
    'type',char(zeros(length(edges_nodes(:,1)),1)),'order',...
    zeros(length(edges_nodes(:,1)),1)); 
Edges.type(:) = 'D';           % all edges are predefined as Dirichlet 
Edges.order(:) = NaN;          % all degrees are predefined as NaN


Loops=struct('nodes',loops_nodes,'edges',loops_edges,...
    'center',zeros(length(loops_nodes(:,1)),2),...
    'area',zeros(length(loops_nodes(:,1)),1),...
    'order',zeros(length(loops_nodes(:,1)),1),...
    'materials',zeros(length(loops_nodes(:,1)),5));

% MATERIAL DATA

load('StructDef','PlaneState','Young','Poisson');
PS = PlaneState;

% Allocate the material characteristics defined in the GUI to all elements
Loops.materials(:,1) = Young;
Loops.materials(:,2) = Poisson;
% ... or overwrite the material characteristics manually
% Loops.materials(5,1) = 1000;      <------- Examples
% Loops.materials(5,2) = 0.2;      <------- Examples


for ii=1:length(Loops.area)
    if (PS == 1) 
        Loops.materials(ii,3) = Loops.materials(ii,1)/...
            (1-(Loops.materials(ii,2))^2);
        Loops.materials(ii,4) = (Loops.materials(ii,2))*...
            ((Loops.materials(ii,1))/(1-(Loops.materials(ii,2))^2));
        Loops.materials(ii,5) = ((1-Loops.materials(ii,2))/2)*...
            (Loops.materials(ii,1)/(1-(Loops.materials(ii,2))^2));
    elseif (PS == 2)
        Loops.materials(ii,3)= (1-Loops.materials(ii,2))...
            *(Loops.materials(ii,1)/((1+Loops.materials(ii,2))...
            *(1-2*Loops.materials(ii,2))));
        Loops.materials(ii,4)= (Loops.materials(ii,2))*...
            (Loops.materials(ii,1)/((1+Loops.materials(ii,2))*...
            (1-2*Loops.materials(ii,2))));
        Loops.materials(ii,5)=((1-2*Loops.materials(ii,2))/2)*...
            (Loops.materials(ii,1)/((1+Loops.materials(ii,2))...
            *(1-2*Loops.materials(ii,2))));
    else
        error('local:consistencyChk',...
            'Ill defined plane state PS = %d. \n',PS);
    end
end

for i=1:length(loops_nodes(:,1))
    [Loops.center(i,:),Loops.area(i)] = ...
        polygonCentroid(Nodes(loops_nodes(i,:),:));
    Loops.area(i)=abs(Loops.area(i)); % for the area to be correct, the 
                             % nodes in loops_nodes must be listed in CW or
                             % CCW order!
end

%% *************************************************
% *************** USER-DEFINED AREA *************** 

% RUN CONTROL DATA
load('StructDef','NumberGaussPoints','EdgesOrder','LoopsOrder');
load('StructRegBC2','edgesDirichlet','edgesNeumann','dataDir','dataNeu');
NGP = NumberGaussPoints;   %Number of Gauss integration points per interval

% EDGE TYPE DATA

% Registration of the Neumann edges, using edgesNeumann vector from the GUI
if exist('edgesNeumann')
    for i=1:length(edgesNeumann)
        Edges.type(edgesNeumann(i),1) = 'N'; 
    end
end
% Manual definition of Neumann edges, if needed
%Edges.type(17:20) = 'N';       <------- Examples       
%Edges.type(25:5:40) = 'N';     <------- Examples

% EDGE REFINEMENT DATA

% Allocate the refinement order defined in the GUI to all Dirichlet 
% boundaries
Edges.order(Edges.type=='D') = EdgesOrder;
% ... or overwrite orders manually, if you need to have different orders
% Edges.order(13) = 3;       <------- Examples
% Edges.order(17) = 5;       <------- Examples

% ELEMENT DATA

% Allocate the refinement order defined in the GUI to all elements
Loops.order(:) = LoopsOrder;
% ... or overwrite orders manually, if you need to have different orders
% for different elements
% Loops.order(1) = 5;       <------- Examples
% Loops.order(2) = 5;       <------- Examples
% Loops.order(3) = 5;       <------- Examples

% LOAD & DISPLACEMENT DATA (evenly spaced points for polynomial interpolation)
% initialization...

BConds=struct('Neumann',{cell(length(edges_nodes),2)},'Dirichlet',...
    {cell(length(edges_nodes),2)});
BConds.Neumann(:) = {NaN};
BConds.Dirichlet(:) = {NaN};

% Dirichlet boundary conditions, as imported from the GUI
if exist('dataDir')
    for i=1:size(dataDir,1)
        BConds.Dirichlet{dataDir{i,1},1}=str2num(dataDir{i,2});
        BConds.Dirichlet{dataDir{i,1},2}=str2num(dataDir{i,3});
    end
end
% overwrite boundary conditions manually
% BConds.Dirichlet(17,1) = {[0 0.25]};  <------- Examples
% BConds.Dirichlet(17,2) = {[0 0.25]};  <------- Examples

% Neumann boundary conditions
if exist('dataNeu')
    for i=1:size(dataNeu,1)
        BConds.Neumann{dataNeu{i,1},1}=str2num(dataNeu{i,2});
        BConds.Neumann{dataNeu{i,1},2}=str2num(dataNeu{i,3});
    end
end
% overwrite boundary conditions manually
% BConds.Neumann(17,1) = {[0 0.25]};  <------- Examples
% BConds.Neumann(17,2) = {[0 0.25]};  <------- Examples

% GRID POINTS FOR THE POST-PROCESSING. 
NoDiv = NGP;

% *************** END OF USER-DEFINED AREA *************** 
% ********************************************************

%% ********* CHECKS **************
for ii=1:length(Edges.type)
    if Edges.type(ii)=='D' && any(~isnan(cell2mat(BConds.Neumann(ii,1)))) ...
            && any(~isnan(cell2mat(BConds.Neumann(ii,2))))
        % if Dirichlet boundary and any of the Neumann BC is not NaN
        error('local:consistencyChk',...
            'Edge %d is defined as Dirichlet, but has Neumann boundary conditions in both directions. \n',...
            ii);
    elseif(Edges.type(ii)=='D' && isnan(Edges.order(ii)))
        % if Dirichlet boundary and the order is NaN
        error('local:consistencyChk',...
            'Edge %d is defined as Dirichlet, but has no order of approximation. \n',...
            ii);
    elseif Edges.type(ii)=='N' &&( any(~isnan(cell2mat(BConds.Dirichlet(ii,1)))) ...
            || any(~isnan(cell2mat(BConds.Dirichlet(ii,2)))))
        % if Neumann boundary and any of the Dirichlet BC is not NaN
        error('local:consistencyChk',...
            'Edge %d is defined as Neumann, but has Dirichlet boundary conditions. \n',...
            ii);
    elseif(Edges.type(ii)=='N' && ~isnan(Edges.order(ii)))
        % if Neumann boundary and the order is not NaN
        error('local:consistencyChk',...
            'Edge %d is defined as Neumann, but has non-zero order of approximation. \n',...
            ii);
    end
end
    
end
