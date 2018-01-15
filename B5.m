function LHS = B5(Edges, Loops, LHS, abscissa, weight)
% B5 sweeps through the edges and calls the functions that generate the
% B5 block of the boundary matrix in the LHS.
%
% B5 is called by MAIN***. It receives as input data the Edges and Loops
% structures, the LHS matrix (that is, the matrix of coefficients of the
% solving system), and the Gauss-Legendre integration parameters abscissa
% and weights. It returns to MAIN*** the LHS matrix with the B5 blocks of
% all elements inserted at the correct positions (as determined in
% ASSIGNPARTS).
%
%
% BIBLIOGRAPHY
% 1. FreeHyTE Page - https://sites.google.com/site/ionutdmoldovan/freehyte
% 2. Moldovan ID, Cismasiu I - FreeHyTE: theoretical bases and developer’s 
% manual, https://drive.google.com/file/d/0BxuR3pKS2hNHTzB5N2Q4cXZKcGc/view
% 3. FreeHyTE Structural HTD User's Manual - 
%    https://drive.google.com/drive/folders/0BxuR3pKS2hNHWWkyb2xmTGhLa1k
% 4. Silva V - Elementos finitos híbridos-Trefftz de deslocamento para 
% problemas de elasticidade plana, MSc Thesis, Universidade Nova de Lisboa,
% 2016 (in Portuguese).
%
%
% B5 computes the internal product between the following bases, expressed
% in a polar (r,th) referential:
% * the rigid body displacement basis U5,
%        U5 = | Ur | = |   0    cos(th)   sin(th) |
%             | Ut |   |   r   -sin(th)   cos(th) |
% * the boundary traction basis Z, defined by Chebyshev polynomials,
%        Z  = cos(m*ArcCos(abscissa))
% where m is the column of the current term in matrix B5.
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
                'insert',Loops.insert(id,5),'dim',Loops.dim(id,5),...
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
                    % block for the current edge and the first rigid body
                    % mode of its left element, in the jj direction.
                    B1ri = sign*B1r_Matrix_i(LocEdge, LocLoop, jj, abscissa, weight);
                    % Inserting the block in the global LHS matrix. The
                    % insertion is made at line LocLoop.insert (the first 
                    % rigid body mode line)& column LocEdge.insert(jj).
                    LHS(LocLoop.insert,...
                        LocEdge.insert(jj):LocEdge.insert(jj)+LocEdge.dim(jj)-1) = -B1ri;
                    % Inserting the conjugate transposed in the global LHS matrix
                    LHS(LocEdge.insert(jj):LocEdge.insert(jj)+LocEdge.dim(jj)-1,...
                        LocLoop.insert) = -B1ri';
 
                    % Launching the function that computes the boundary
                    % block for the current edge and the second rigid body
                    % mode of its left element, in the jj direction.
                    B2ri = sign*B2r_Matrix_i(LocEdge, LocLoop, jj, abscissa, weight);
                    % Inserting the block in the global LHS matrix. The
                    % insertion is made at line LocLoop.insert+1 (the 2nd 
                    % rigid body mode line)& column LocEdge.insert(jj).
                    LHS(LocLoop.insert +1 ,...
                        LocEdge.insert(jj):LocEdge.insert(jj)+LocEdge.dim(jj)-1) = -B2ri;
                    % Inserting the conjugate transposed in the global LHS matrix
                    LHS(LocEdge.insert(jj):LocEdge.insert(jj)+LocEdge.dim(jj)-1,...
                        LocLoop.insert+1 ) = -B2ri';
                    
                    % Launching the function that computes the boundary
                    % block for the current edge and the third rigid body
                    % mode of its left element, in the jj direction.
                    B3ri = sign*B3r_Matrix_i(LocEdge, LocLoop, jj, abscissa, weight);
                    % Inserting the block in the global LHS matrix. The
                    % insertion is made at line LocLoop.insert+2 (the 3rd 
                    % rigid body mode line)& column LocEdge.insert(jj).
                    LHS(LocLoop.insert+2,...
                        LocEdge.insert(jj):LocEdge.insert(jj)+LocEdge.dim(jj)-1) = -B3ri;
                    % Inserting the conjugate transposed in the global LHS matrix
                    LHS(LocEdge.insert(jj):LocEdge.insert(jj)+LocEdge.dim(jj)-1,...
                        LocLoop.insert+2) = -B3ri';
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
                'insert',Loops.insert(id,5),'dim',Loops.dim(id,5));
            
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
                    % block for the current edge and the first rigid body
                    % mode of its right element, in the jj direction.
                    B1ri = sign*B1r_Matrix_i(LocEdge, LocLoop, jj, abscissa, weight);
                    % Inserting the block in the global LHS matrix. The
                    % insertion is made at line LocLoop.insert (the first 
                    % rigid body mode line)& column LocEdge.insert(jj).
                    LHS(LocLoop.insert,...
                        LocEdge.insert(jj):LocEdge.insert(jj)+LocEdge.dim(jj)-1) = -B1ri;
                    % Inserting the conjugate transposed in the global LHS matrix
                    LHS(LocEdge.insert(jj):LocEdge.insert(jj)+LocEdge.dim(jj)-1,...
                        LocLoop.insert) = -B1ri';

                    % Launching the function that computes the boundary
                    % block for the current edge and the second rigid body
                    % mode of its right element, in the jj direction.
                    B2ri = sign*B2r_Matrix_i(LocEdge, LocLoop, jj, abscissa, weight);
                    % Inserting the block in the global LHS matrix. The
                    % insertion is made at line LocLoop.insert+1 (the 2nd
                    % rigid body mode line)& column LocEdge.insert(jj).
                    LHS(LocLoop.insert+1,...
                        LocEdge.insert(jj):LocEdge.insert(jj)+LocEdge.dim(jj)-1) = -B2ri;
                    % Inserting the conjugate transposed in the global LHS matrix
                    LHS(LocEdge.insert(jj):LocEdge.insert(jj)+LocEdge.dim(jj)-1,...
                        LocLoop.insert+1) = -B2ri';
                    
                    % Launching the function that computes the boundary
                    % block for the current edge and the third rigid body
                    % mode of its right element, in the jj direction.
                    B3ri = sign*B3r_Matrix_i(LocEdge, LocLoop, jj, abscissa, weight);
                    % Inserting the block in the global LHS matrix. The
                    % insertion is made at line LocLoop.insert+2 (the 3rd
                    % rigid body mode line)& column LocEdge.insert(jj).
                    LHS(LocLoop.insert+2,...
                        LocEdge.insert(jj):LocEdge.insert(jj)+LocEdge.dim(jj)-1) = -B3ri;
                    % Inserting the conjugate transposed in the global LHS matrix
                    LHS(LocEdge.insert(jj):LocEdge.insert(jj)+LocEdge.dim(jj)-1,...
                        LocLoop.insert+2) = -B3ri';
                end
            end
        end
    end
end
end



function B1ri = B1r_Matrix_i(LocEdge, LocLoop, jj, abscissa, weight)
% B1R_MATRIX_I local function computes the B5 boundary block for the LocEdge
% edge and the first rigid body mode of the LocLoop element, in the jj 
% direction. The side is mapped to a [-1,1] interval to perform the 
% integration.

%% Initialization 
% Initialization of the B5 (R1) block
B1ri = zeros(LocLoop.dim,LocEdge.dim(jj));

% a single rigid body mode; m+1 is the current column
n = 1; 
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
Ur = 0;
Ut = R;

% Computing the normal or tangential projection of the displacement basis,
% according to the current direction, jj
if jj == 1  % normal direction
    U = NR.*Ur + NT.*Ut;
else     % tangential direction
    U = -NT.*Ur + NR.*Ut;
end

% Computing the traction basis - Chebyshev polynomials mapped on the
% abscissa [-1,1] interval are used.
Z = cos(bsxfun(@times,M,acos(A)));

%% Computing the boundary integral
% The integral is the internal product between the displacement basis U
% (projected in the direction jj) and the traction basis Z
B1ri3D = U.*Z;
w3D(1,1,:) = weight;
B1ri = L/2 * sum(bsxfun(@times,B1ri3D,w3D),3); % computes the integral

end



function B2ri = B2r_Matrix_i(LocEdge, LocLoop, jj, abscissa, weight)
% B2R_MATRIX_I local function computes the B5 boundary block for the LocEdge
% edge and the second rigid body mode of the LocLoop element, in the jj 
% direction. The side is mapped to a [-1,1] interval to perform the 
% integration.

%% Initialization 
% Initialization of the B5 (R2) block
B2ri = zeros(LocLoop.dim,LocEdge.dim(jj));

% a single rigid body mode; m+1 is the current column
n = 1;
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
Ur = cos(T);
Ut = -sin(T);
% Computing the normal or tangential projection of the displacement basis,
% according to the current direction, jj
if jj == 1  % normal direction
    U = NR.*Ur + NT.*Ut;
else     % tangential direction
    U = -NT.*Ur + NR.*Ut;
end

% Computing the traction basis - Chebyshev polynomials mapped on the
% abscissa [-1,1] interval are used.
Z = cos(bsxfun(@times,M,acos(A)));

%% Computing the boundary integral
% The integral is the internal product between the displacement basis U
% (projected in the direction jj) and the traction basis Z
B2ri3D = U.*Z;
w3D(1,1,:) = weight;
B2ri = L/2 * sum(bsxfun(@times,B2ri3D,w3D),3); % computes the integral

end



function B3ri = B3r_Matrix_i(LocEdge, LocLoop, jj, abscissa, weight)
% B3R_MATRIX_I local function computes the B5 boundary block for the LocEdge
% edge and the third rigid body mode of the LocLoop element, in the jj 
% direction. The side is mapped to a [-1,1] interval to perform the 
% integration.

%% Initialization 
% Initialization of the B5 (R3) block
B3ri = zeros(LocLoop.dim,LocEdge.dim(jj));

% a single rigid body mode; m+1 is the current column
n = 1;
m = 0:LocEdge.order(1);                      

%% Generating the geometric data
% The following code transforms the abscissa coordinates, expressed in
% the [-1,1] referential, to the polar coordinates required to compute
% the values of the basis functions. The components of the outward
% normal to the boundary in the radial and tangential directions are
% also calculated. They are required to compute the normal and
% tangential components of the displacement basis.

% Computing the length of the current edge
L = sqrt(LocEdge.parametric(3)^2 + LocEdge.parametric(4)^2); 

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
Ur = sin(T);
Ut = cos(T);

% Computing the normal or tangential projection of the displacement basis,
% according to the current direction, jj
if jj == 1  % normal direction
    U = NR.*Ur + NT.*Ut;
else     % tangential direction
    U = -NT.*Ur + NR.*Ut;
end

% Computing the traction basis - Chebyshev polynomials mapped on the
% abscissa [-1,1] interval are used.
Z = cos(bsxfun(@times,M,acos(A)));

%% Computing the boundary integral
% The integral is the internal product between the displacement basis U
% (projected in the direction jj) and the traction basis Z
B3ri3D = U.*Z;
w3D(1,1,:) = weight;
B3ri = L/2 * sum(bsxfun(@times,B3ri3D,w3D),3); % computes the integral

end