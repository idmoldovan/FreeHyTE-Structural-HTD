function LHS = B3(Edges, Loops, LHS, abscissa, weight)
% B3 sweeps through the edges and calls the functions that generate the
% B3 block of the boundary matrix in the LHS.
%
% B3 is called by MAIN***. It receives as input data the Edges and Loops
% structures, the LHS matrix (that is, the matrix of coefficients of the
% solving system), and the Gauss-Legendre integration parameters abscissa
% and weights. It returns to MAIN*** the LHS matrix with the B3 blocks of
% all elements inserted at the correct positions (as determined in
% ASSIGNPARTS).
%
%
% BIBLIOGRAPHY
% 1. FreeHyTE Page - https://sites.google.com/site/ionutdmoldovan/freehyte
% 2. Moldovan ID, Cismasiu I - FreeHyTE: theoretical bases and developer�s 
% manual, https://drive.google.com/file/d/0BxuR3pKS2hNHTzB3N2Q4cXZKcGc/view
% 3. FreeHyTE Structural HTD User's Manual - 
%    https://drive.google.com/drive/folders/0BxuR3pKS2hNHWWkyb2xmTGhLa1k
% 4. Silva V - Elementos finitos h�bridos-Trefftz de deslocamento para 
% problemas de elasticidade plana, MSc Thesis, Universidade Nova de Lisboa,
% 2016 (in Portuguese).
%
%
% B3 computes the internal product between the following bases, expressed
% in a polar (r,th) referential:
% * the displacement basis U3, generated by odd biharmonic potentials,
%        U3 = | Ur | = r^n * | [k3*(n-3) + k2*(n-1)] * sin[(n-1)*th] |
%             | Ut |         | [k3*(n+3) + k2*(n+1)] * cos[(n-1)*th] |
% * the boundary traction basis Z, defined by Chebyshev polynomials,
%        Z  = cos(m*ArcCos(abscissa))
% where n and m are the line and column of the current term in matrix B3,
% and k2 and k3 are stiffness coefficients (see % INPUTPROC***).
%
% Consistent with the direction of the boundary traction approximation, the 
% displacement basis is projected onto the boundary normal or tangential
% direction,
%        Un = nr*Ur + nt*Ut
%        Ut = -nt*Ur + nr*Ut
% where nr and nt are the radial and tangential components of the outward 
% unit normal to the current boundary.
%
%
% Further details on the structure os the solving system are presented in 
% reference [2] (Section 6.2). The transformation of referentials is
% covered in reference [4] (Appendix C). 

%% Sweeping the edges
for ii=1:length(Edges.type)
    % Boundary blocks are constructed for Dirichlet (and interior)
    % boundaries only.
    if strcmpi(Edges.type(ii),'D')
        % LocEdge is a local structure where the features of the current
        % edge which are directly useful for the calculation of the
        % boundary block are stored.
        LocEdge =struct('nini',Edges.nini(ii),'nfin',Edges.nfin(ii),...
            'parametric',Edges.parametric(ii,:), 'lleft',Edges.lleft(ii),...
            'lright',Edges.lright(ii), 'order',Edges.order(ii),...
            'insert',Edges.insert(ii,:), 'dim',Edges.dim(ii,:));
        
        % Generating the boundary block corresponding to the current edge
        % and its LEFT finite element
        if LocEdge.lleft
            id = LocEdge.lleft;
            sign = 1.;
            % LocLoop is a local structure where the features of the left
            % element which are directly useful for the calculation of the
            % boundary block are stored.
            LocLoop = struct('id',id,'edges',Loops.edges(id,:),... 
                'center',Loops.center(id,:),'order',Loops.order(id,1),...
                'insert',Loops.insert(id,3),'dim',Loops.dim(id,3),...
                'materials',Loops.materials(id,:));
            
            % The boundary blocks are generated separately in the normal
            % and tangential directions, if approximations are present. 
            % jj = 1 corresponds to the normal direction. 
            % jj = 2 denotes the tangential direction.
            for jj = 1:2
                % Verifying if an approximation is present in the jj
                % direction. If no approximations are present,
                % LocEdge.dim(jj) is null (see ASSIGNPARTS).
                if LocEdge.dim(jj)
                    % Launching the function that computes the boundary
                    % block for the current edg and its left element, in 
                    % the jj direction.
                    B3i = sign*B3_Matrix_i(LocEdge, LocLoop, jj, abscissa, weight);
                    % Inserting the block in the global LHS matrix. The 
                    % insertion is made at line LocLoop.insert & column 
                    % LocEdge.insert(jj).
                    LHS(LocLoop.insert:LocLoop.insert+LocLoop.dim-1,...
                        LocEdge.insert(jj):LocEdge.insert(jj)+LocEdge.dim(jj)-1) = -B3i;
                    % Inserting the conjugate transposed in the global LHS matrix
                    LHS(LocEdge.insert(jj):LocEdge.insert(jj)+LocEdge.dim(jj)-1,...
                        LocLoop.insert:LocLoop.insert+LocLoop.dim-1) = -B3i';
                end
            end
        else                       % there should always be a left element
            error('local:consistencyChk',...
                'No left loop for edge %d. \n', ii);
        end
        
        % Generating the boundary block corresponding to the current edge
        % and its RIGHT finite element 
        if LocEdge.lright
            id = LocEdge.lright;
            sign = +1.;
            % LocLoop is a local structure where the features of the right
            % element which are directly useful for the calculation of the
            % boundary block are stored.
            LocLoop = struct('id',id,'edges',Loops.edges(id,:),...
                'center',Loops.center(id,:),'order',Loops.order(id,1),...
                'insert',Loops.insert(id,3),'dim',Loops.dim(id,3),...
                'materials',Loops.materials(id,:));
            % The boundary blocks are generated separately in the normal
            % and tangential directions, if approximations are present.
            % jj = 1 corresponds to the normal direction.
            % jj = 2 denotes the tangential direction.
            for jj = 1:2
                % Verifying if an approximation is present in the jj
                % direction. If no approximations are present,
                % LocEdge.dim(jj) is null (see ASSIGNPARTS).
                if LocEdge.dim(jj)
                    % Launching the function that computes the boundary
                    % block for the current edge and its right element, in 
                    % the jj direction.
                    B3i = sign*B3_Matrix_i(LocEdge, LocLoop, jj, abscissa, weight);
                    % Inserting the block in the global LHS matrix. The 
                    % insertion is made at line LocLoop.insert & column 
                    % LocEdge.insert(jj).
                    LHS(LocLoop.insert:LocLoop.insert+LocLoop.dim-1,...
                        LocEdge.insert(jj):LocEdge.insert(jj)+LocEdge.dim(jj)-1) = -B3i;
                    % Inserting the conjugate transposed in the global LHS matrix
                    LHS(LocEdge.insert(jj):LocEdge.insert(jj)+LocEdge.dim(jj)-1,...
                        LocLoop.insert:LocLoop.insert+LocLoop.dim-1) = -B3i';
                end
            end
        end
    end
end
end

function B3i = B3_Matrix_i(LocEdge, LocLoop, jj, abscissa, weight)
% B3_MATRIX_I local function computes the B3 boundary block for the LocEdge
% edge and LocLoop element, in the jj direction. The side is mapped to a 
% [-1,1] interval to perform the integration.

%% Initialization 
% Initialization of the B3 block
B3i = zeros(LocLoop.dim,LocEdge.dim(jj));

% n-1 denotes the current line (the term of the displacement basis is a
% rigid body term for n = 1); m+1 is the current column
n = 2:LocLoop.order(1);
m = 0:LocEdge.order(1);                      

%% Generating the geometric data
% The following code transforms the abscissa coordinates, expressed in
% the [-1,1] referential, to the polar coordinates required to compute
% the values of the basis functions. The components of the outward
% normal to the boundary in the radial and tangential directions are
% also calculated. They are required to compute the normal and
% tangential components of the displacement basis.

% Computing the length of the current edge
L = sqrt(LocEdge.parametric(3)^2 + LocEdge.parametric(4)^2); % length

% Constructing the 3D matrices containing the n x m x abscissa
% integration grid
[N,M,A] = ndgrid(n,m,abscissa);

% Transforming the edge abscissa into local coordinates. The local
% referential is centered in the barycenter of the element, its axes
% aligned with the Cartesian axes of the global referential.
loc_x = LocEdge.parametric(1) - LocLoop.center(1) + 0.5 *...
    (A + 1) * LocEdge.parametric(3);  
loc_y = LocEdge.parametric(2) - LocLoop.center(2) + 0.5 *...
    (A + 1) * LocEdge.parametric(4);

% Transforming the local Cartesian coordinates into polar.
R = sqrt(loc_x.^2 + loc_y.^2);  
T = atan2(loc_y, loc_x);

% Computing the components of the outward normal in the Cartesian
% directions.
nx = LocEdge.parametric(4) / L;
ny = -1* LocEdge.parametric(3) / L;
if LocEdge.lright==LocLoop.id  % if the element is on the right,
    nx = -nx;                  % change the sign of the normal
    ny = -ny;
end

% Computing the components of the outward normal in the polar directions.
NR = nx * cos(T) + ny * sin(T);   
NT = -1*nx * sin(T) + ny * cos(T);

%% Computing the basis functions for all integration points
% Polar components of the displacement basis
Ur = R.^N.*sin((N-1).*T).*((LocLoop.materials(5)*(N-3))+...
    (LocLoop.materials(4)*(N-1))) ;
Ut = R.^N.*cos((N-1).*T).*((LocLoop.materials(5)*(N+3))+...
    (LocLoop.materials(4)*(N+1))) ;

% Computing the normal or tangential projection of the displacement basis,
% according to the current direction, jj
if jj == 1  % normal direction
    U = NR.*Ur + NT.*Ut;
else    % tangential direction
    U = -NT.*Ur + NR.*Ut;
end

% Computing the traction basis - Chebyshev polynomials mapped on the
% abscissa [-1,1] interval are used.
Z = cos(bsxfun(@times,M,acos(A)));

%% Computing the boundary integral
% The integral is the internal product between the displacement basis U
% (projected in the direction jj) and the traction basis Z
B3i3D = U.*Z;
w3D(1,1,:) = weight;
B3i = L/2 * sum(bsxfun(@times,B3i3D,w3D),3); % computes the integral

end