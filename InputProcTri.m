function [NGP, Nodes, Edges, Loops, BConds, NoDiv] =  InputProcTri
% INPUTPROCTRI is the input processing function for non-regular meshes.
%
% INPUTPROCTRI is called by MAINTRI. It reads the input data inserted in 
% the GUIs and organizes it in the data structures Edges, Loops and BConds.
%
%
% BIBLIOGRAPHY
% 1. FreeHyTE Page - https://sites.google.com/site/ionutdmoldovan/freehyte
% 2. Moldovan ID, Cismasiu I - FreeHyTE: theoretical bases and developer’s 
% manual, https://drive.google.com/file/d/0BxuR3pKS2hNHTzB2N2Q4cXZKcGc/view
% 3. FreeHyTE Structural HTD User's Manual - 
%    https://drive.google.com/drive/folders/0BxuR3pKS2hNHWWkyb2xmTGhLa1k
% 4. Silva V - Elementos finitos híbridos-Trefftz de deslocamento para 
% problemas de elasticidade plana, MSc Thesis, Universidade Nova de Lisboa,
% 2016 (in Portuguese).
%
%
% INPUT DATA (read from the *.mat files created by the GUIs):
% * from STRUCTDEF: PlaneState (is equal to 1 for plane stress and 2 for 
% plane strain), Young (Young's modulus), Poisson (Poisson coefficient),
% NumberGaussPoints (the nuber of Gauss-Legendre points for the side
% integration), EdgesOrder and LoopsOrder (the orders of the approximation
% bases on the essential edges of the mesh and in the finite elements);
% * from STRUCTTRIBC1: nodes, edges_nodes, edges_loops, loops_nodes and 
% loops_edges. These variables contain geometrical and topological 
% information regarding the mesh. For a full description of these
% variables, please consult Section 5.3 of reference [2];
% * from STRUCTTRIBC2: edgesDirichlet, edgesNeumann, dataDir, dataNeu;
% * edgesDirichlet (edgesNeumann): list with the Dirichlet (Neumann) edges;
% * dataDir (dataNeu): cell array with three columns and as many lines as
% Dirichlet (Neumann) edges. For each edge, it stores the edge's index, 
% and the Dirichlet (Neumann) boundary conditions in the normal and 
% tangential directions. The boundary conditions are stored as strings. For
% further details regarding the definition of the boundary conditions,
% please refer to Section 4.6 of reference [3].
%
%
% OUTPUT DATA to function MAINTRI:
% * NGP is the number of Gauss points for the line integration;
% * Nodes is a (NNODE x 2) matrix, where NNODE is the number of nodes in 
% mesh. It stores the coordinates of each node;
% * Edges, Loops and BConds are data structures storing information
% on the edges, finite elements (loops) and boundary conditions,
% respectively. They are documented in Section 5.3 of reference [2];
% * NoDiv^2 is the number of points for plotting the colormaps of the
% solution.

%% MESH DATA
% Creation of the Nodes, Edges and Loops data structures

% Loading mesh information from STRUCTDEF GUI
load('StructTriBC1','nodes','edges_nodes','edges_loops','loops_nodes',...
     'loops_edges');

% Definition of the mesh-related data Nodes, and data structures Edges and
% Loops. For a full description of these data structure, please refer to
% Section 5.3 of reference [2].
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

% Computation of the barycenter and area of each finite element. It uses
% the POLYGONCENTROID function by David Legland.
for i=1:length(loops_nodes(:,1))
    [Loops.center(i,:),Loops.area(i)] = ...
        polygonCentroid(Nodes(loops_nodes(i,:),:));
    Loops.area(i)=abs(Loops.area(i));
end

%% MATERIAL DATA
% Loads the material data and stores them in the Loop structure. The
% material data is assumed uniform for all elements. However, users may
% overwrite the data loaded from the GUI to define elements with distinct
% material properties. An example is given below.

% Loading the type of plane problem and the material data.
load('StructDef','PlaneState','Young','Poisson');
PS = PlaneState; % 1 = plane stress; 2 = plane strain.

% Allocate the material characteristics defined in the GUI to all elements
Loops.materials(:,1) = Young;
Loops.materials(:,2) = Poisson;
% ... or overwrite the material characteristics manually
% Loops.materials(5,1) = 1000;      <------- Examples
% Loops.materials(5,2) = 0.2;      <------- Examples

% Computing the three material stiffness coefficients identified below in
% the generic Hooke's law for plane elasticity:
%
%  |sxx|  =  |k1  k2  0 | * |exx|
%  |syy|     |k2  k1  0 |   |eyy|
%  |sxy|     |0   0   k3|   |exy|
%
% The expressions of the stiffness coefficients depend on the type of plane
% state.
for ii=1:length(Loops.area)
    if (PS == 1) % Plane stress
        % Computing k1
        Loops.materials(ii,3) = Loops.materials(ii,1)/...
            (1-(Loops.materials(ii,2))^2);
        % Computing k2
        Loops.materials(ii,4) = (Loops.materials(ii,2))*...
            ((Loops.materials(ii,1))/(1-(Loops.materials(ii,2))^2));
        % Computing k3
        Loops.materials(ii,5) = ((1-Loops.materials(ii,2))/2)*...
            (Loops.materials(ii,1)/(1-(Loops.materials(ii,2))^2));
    else         % Plane strain
        % Computing k1 
        Loops.materials(ii,3)= (1-Loops.materials(ii,2))...
            *(Loops.materials(ii,1)/((1+Loops.materials(ii,2))...
            *(1-2*Loops.materials(ii,2))));
        % Computing k2 
        Loops.materials(ii,4)= (Loops.materials(ii,2))*...
            (Loops.materials(ii,1)/((1+Loops.materials(ii,2))*...
            (1-2*Loops.materials(ii,2))));
        % Computing k3 
        Loops.materials(ii,5)=((1-2*Loops.materials(ii,2))/2)*...
            (Loops.materials(ii,1)/((1+Loops.materials(ii,2))...
            *(1-2*Loops.materials(ii,2))));
    end
end
    

%% RUN CONTROL DATA
% Loading algorithmic, refinement and edge data
load('StructDef','NumberGaussPoints','EdgesOrder','LoopsOrder');
load('StructTriBC2','edgesDirichlet','edgesNeumann','dataDir','dataNeu');
NGP = NumberGaussPoints;   %Number of Gauss integration points per interval

%% EDGE TYPE DATA
% Registration of the Neumann edges, using edgesNeumann vector from the
% GUI. It is recalled that all edges were predefined as Dirichlet.
% Users may overwrite the data loaded from the GUI to change the boundary
% types. An example is given below.
if exist('edgesNeumann')
    for i=1:length(edgesNeumann)
        Edges.type(edgesNeumann(i),1) = 'N'; 
    end
end
% Manual definition of Neumann edges, if needed
%Edges.type(17:20) = 'N';       <------- Examples       
%Edges.type(25:5:40) = 'N';     <------- Examples

%% EDGE REFINEMENT DATA
% Allocates the refinement order defined in the GUI to all Dirichlet 
% boundaries
Edges.order(Edges.type=='D') = EdgesOrder;
% ... or overwrite orders manually, if you need to have different orders
% Edges.order(13) = 3;       <------- Examples
% Edges.order(17) = 5;       <------- Examples

%% ELEMENT DATA
% Allocates the refinement order defined in the GUI to all elements
Loops.order(:) = LoopsOrder;
% ... or overwrite orders manually, if you need to have different orders
% for different elements
% Loops.order(1) = 5;       <------- Examples
% Loops.order(2) = 5;       <------- Examples
% Loops.order(3) = 5;       <------- Examples

%% LOAD & DISPLACEMENT DATA 
% Boundary conditions can be described by polynomials of any order. The 
% definition of a boundary condition is made by specifying its values in 
% as many equally spaced points along the boundary as needed to define its
% polynomial variation. For further details regarding the definition of the 
% boundary conditions, please refer to Section 4.6 of reference [3].
%
% BConds data structure collects information regarding the boundary 
% conditions. Its members are cell arrays with as many lines as the 
% external boundaries of the structure. The values of the forces (or
% displacements) enforced on the boundaries are stored in the Neumann (or
% Dirichlet) fields of the structure. NaN is stored in the Dirichlet field
% of a Neumann boundary, and vice-versa.

% Initialization of the BConds structure
BConds=struct('Neumann',{cell(length(edges_nodes),2)},'Dirichlet',...
    {cell(length(edges_nodes),2)});
BConds.Neumann(:) = {NaN};
BConds.Dirichlet(:) = {NaN};

% Dirichlet boundary conditions are imported from the GUI and stored in the
% Dirichlet field of the structure, in the normal and tangential
% directions. Users may overwrite the data loaded from the GUI to change 
% the boundary conditions. An example is given below.
if exist('dataDir')
    for i=1:size(dataDir,1)
        BConds.Dirichlet{dataDir{i,1},1}=str2num(dataDir{i,2});
        BConds.Dirichlet{dataDir{i,1},2}=str2num(dataDir{i,3});
    end
end
% overwrite boundary conditions manually
% BConds.Dirichlet(17,1) = {[0 0.25]};  <------- Examples
% BConds.Dirichlet(17,2) = {[0 0.25]};  <------- Examples

% Neumann boundary conditions are imported from the GUI and stored in the
% Neumann field of the BConds structure, in the normal and tangential
% directions. Users may overwrite the data loaded from the GUI to change 
% the boundary conditions. An example is given below.
if exist('dataNeu')
    for i=1:size(dataNeu,1)
        BConds.Neumann{dataNeu{i,1},1}=str2num(dataNeu{i,2});
        BConds.Neumann{dataNeu{i,1},2}=str2num(dataNeu{i,3});
    end
end
% overwrite boundary conditions manually
% BConds.Neumann(17,1) = {[0 0.25]};  <------- Examples
% BConds.Neumann(17,2) = {[0 0.25]};  <------- Examples

% NoDiv is the number of points for plotting the colormaps of the
% solution, in each Cartesian direction. It is predefined here as the
% number of Gauss points used for the side integration, but it may be
% arbitrarily defined to some other value.
NoDiv = NGP;


%% CHECKS
% The consistency of the boundary definition is checked here.
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
