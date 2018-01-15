function RHS = U(Edges, BConds, RHS, abscissa, weight)
% U sweeps through the exterior Dirichlet edges and calls the functions 
% generates the U blocks of the free vector (RHS) of the solving system. 
% Neumann edges are not associated with U blocks in the free vector. For 
% interior edges, the U blocks are filled with zeros.
%
% U is called by MAIN***. It receives as input data the Edges, Loops and
% BConds structures, the RHS vector (that is, the free vector of the
% solving system), and the Gauss-Legendre integration parameters abscissa
% and weights. It returns to MAIN*** the RHS vector with the U blocks of
% all exterior Dirichlet edges inserted at the correct positions (as 
% determined in ASSIGNPARTS).
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
% U computes the internal product between the boundary traction basis Z of 
% the Dirichlet edge and the applied displacements on the same edge. 
% * the boundary traction basis Z is defined by Chebyshev polynomials,
%        Z  = cos(m*ArcCos(abscissa))
% where m is the line of the current entry in vector U. 
% * the boundary displacement function is defined in the GUI by its values 
% in an arbitrary number of equally-spaced points along the boundary (see
% Section 4.6 of reference [3]). A polynomial interpolation is performed
% between these values to obtain the analytic expression of the applied
% displacements. The degree of the polynomial is equal to the number of
% displacement values, minus one.
%
% Further details on the structure of the solving system are presented in
% reference [2] (Section 6.2).

%% Sweeping the edges and selecting the exterior Dirichlet boundaries
for ii=1:length(Edges.type)
    
    % Exterior Dirichlet boundaries have no right element
    if (strcmpi(Edges.type(ii),'D') && Edges.lright(ii) == 0 )
        % LocEdge is a local structure where the features of the current
        % edge which are directly useful for the calculation of the
        % U block are stored.
        LocEdge =struct('id',ii,'nini',Edges.nini(ii),'nfin',Edges.nfin(ii),...
            'parametric',Edges.parametric(ii,:), 'lleft',Edges.lleft(ii),...
            'lright',Edges.lright(ii), 'order',Edges.order(ii),...
            'insert',Edges.insert(ii,:), 'dim',Edges.dim(ii,:));
        
        % On an exterior Dirichlet boundary, enforced displacements may be
        % applied in one or two directions of the normal-tangential
        % boundary referential. The definition of the referentials is
        % covered in reference [2] (Appendix A). U blocks are generated
        % in the directions where enforced displacements are applied.
        % jj = 1 corresponds to the normal direction.
        % jj = 2 denotes the tangential direction.
        for jj = 1:2
            % Verifying if an approximation is present in the jj
            % direction. If no approximations are present,LocEdge.dim(jj) 
            % is null (see ASSIGNPARTS).
            if LocEdge.dim(jj)
                
                % Computing the U vector of edge ii. Function U_VECTOR_I is
                % a local function (see below).
                Ui = U_Vector_i(LocEdge, BConds, jj, abscissa, weight);
                
                % Inserting the U vector in the global RHS vector. 
                % The insertion is made at line Edges.insert(ii,jj).
                RHS(LocEdge.insert(jj):LocEdge.insert(jj)+...
                    LocEdge.dim(jj)-1) = -Ui;
            end
        end
    end
end
end

function Ui = U_Vector_i(LocEdge, BConds, jj, abscissa, weight)
% U_VECTOR_I local function computes the U vector of the LocEdge exterior
% Dirichlet boundary. The edge is mapped to a [-1,1] interval to perform 
% the integration.

%% Initialization 
% Initialization of the U block
Ui = zeros(LocEdge.dim(jj),1);

% m+1 denotes the current line of the block
m = 0:LocEdge.order;

%% Generating the geometric data
% Computing the length of the current edge
L = sqrt(LocEdge.parametric(3)^2 + LocEdge.parametric(4)^2); 

%% Computing the traction basis at the integration points
Z = cos(bsxfun(@times,m,acos(abscissa)));

% Computing the values of the enforced displacements at the integration
% points (abscissas)

% obtaining the equally spaced points on [-1,1] interval where the
% displacements are defined and stored in BConds.Dirichlet
a = linspace(-1,1,length(BConds.Dirichlet{LocEdge.id,jj}));
% obtaining the polynomial that interpolates the values in BConds.Dirichlet
if (isnan(BConds.Dirichlet{LocEdge.id,jj})) % just a NaN check
    error('local:consistencyChk',...
        'No Dirichlet boundary conditions are defined on edge %d. \n',...
        LocEdge.id(jj));
else
    pol = polyfit(a,BConds.Dirichlet{LocEdge.id,jj},...
        length(BConds.Dirichlet{LocEdge.id,jj})-1);
end

% computing the values of the interpolation polynomials at the abscissas
q = polyval(pol,abscissa);

%% Computing the integral on the side
% The integral is the internal product between the traction basis and the 
% applied displacements in the jj direction
Ui2D = bsxfun(@times, Z, q);

% Performing the side integration 
Ui = L/2 * sum(bsxfun(@times,Ui2D,weight),1); 
Ui = Ui.';

end