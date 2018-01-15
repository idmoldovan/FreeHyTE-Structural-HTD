function LHS = K14(Edges, Loops,  LHS, abscissa, weight)
% K14 sweeps through the elements and calls the functions that generate the
% K14 block of the stiffness matrix in the LHS.
%
% K14 is called by MAIN***. It receives as input data the Edges and Loops
% structures, the LHS matrix (that is, the matrix of coefficients of the
% solving system), and the Gauss-Legendre integration parameters abscissa
% and weights. It returns to MAIN*** the LHS matrix with the K14 blocks of
% all elements inserted at the correct positions (as determined in
% ASSIGNPARTS).
%
%
% BIBLIOGRAPHY
% 1. FreeHyTE Page - https://sites.google.com/site/ionutdmoldovan/freehyte
% 2. Moldovan ID, Cismasiu I - FreeHyTE: theoretical bases and developer�s 
% manual, https://drive.google.com/file/d/0BxuR3pKS2hNHTzB2N2Q4cXZKcGc/view
% 3. FreeHyTE Structural HTD User's Manual - 
%    https://drive.google.com/drive/folders/0BxuR3pKS2hNHWWkyb2xmTGhLa1k
% 4. Silva V - Elementos finitos h�bridos-Trefftz de deslocamento para 
% problemas de elasticidade plana, MSc Thesis, Universidade Nova de Lisboa,
% 2016 (in Portuguese).
%
%
% K14 computes the internal product between the following bases, expressed
% in a polar (r,th) referential:
% * the displacement basis U1, generated by odd harmonic potentials,
%        U1 = | Ur | = r^n * | sin[(n+1)*th] |
%             | Ut |         | cos[(n+1)*th] |
% * the boundary traction basis N * S4, where S4 is the stress basis
% generated by even biharmonic potentials,
%      S4 = | Sr | = 2*k3(k2+k3)* m * r^(m-1) * | (m-3) * cos[(m-1)*th] |
%           | St |                              |-(m+1) * cos[(m-1)*th] |
%           | Srt|                              |-(m-1) * sin[(m-1)*th] |
%   and N is the directory cosine matrix,
%         N = | nr  0  nt |
%             | 0  nt  nr |
% In the above expressions, n and m are the line and column of the current 
% term in matrix K14, k2 and k3 are stiffness coefficients (see 
% INPUTPROC***), and nr and nt are the radial and tangential components of
% the outward unit normal to the current boundary.
%
% As typical of the Trefftz method, the displacement and stress bases solve
% exactly the equilibrium, compatibility and elasticity equations in the
% domain of each finite element. As a direct consequence, the stiffness
% matrix can be calculated using boundary integrals only.
%
% The derivation of the bases from the solution of the (free-field) Navier
% equation is presented in reference [4] (Section 3.3).

%% Sweeping the elements

for ii=1:length(Loops.area)
    
    % LocLoop is a local structure where the features of the current
    % element which are directly useful for the calculation of the
    % stiffness block are stored.
    LocLoop = struct('id',ii,'edges',Loops.edges(ii,:), 'center',...
        Loops.center(ii,:),'order',Loops.order(ii,:),... 
        'insert',Loops.insert(ii,:),'dim',Loops.dim(ii,:),...
        'materials',Loops.materials(ii,:));
    
    % Computing the K14 matrix of element ii. Function K14_MATRIX_I is a
    % local function (see below).
    K14i = K14_Matrix_i(LocLoop, Edges, abscissa, weight);
    
    % Inserting the matrix in the global LHS matrix. The insertion is made
    % at line Loops.insert(ii,1) and column Loops.insert(ii,4).
    LHS(LocLoop.insert(1):LocLoop.insert(1)+LocLoop.dim(1)-1,...
        LocLoop.insert(4):LocLoop.insert(4)+LocLoop.dim(4)-1) = K14i;
    % Inserting the transposed of the matrix in the global LHS matrix
    LHS(LocLoop.insert(4):LocLoop.insert(4)+LocLoop.dim(4)-1,...
        LocLoop.insert(1):LocLoop.insert(1)+LocLoop.dim(1)-1) = K14i';
    
end

end

function K14i = K14_Matrix_i(LocLoop, Edges, abscissa, weight)
% K14_MATRIX_I local function computes the K14 stiffness block of the 
% LocLoop element. The sides are mapped to a [-1,1] interval to perform the 
% integrations.

%% Initialization 
% Initialization of the K14 block
K14i = zeros([LocLoop.dim(1),LocLoop.dim(4)]);

% n denotes the current line; m is the current column
n = 1:LocLoop.order(1);
m = 1:LocLoop.order(1);

% Sweeping the edges for contour integration
for jj = 1:length(LocLoop.edges)  
    
    % identification of the jj-th edge of the loop
    id = LocLoop.edges(jj);  
    
    % LocEdge is a local structure where the features of the current
    % edge which are directly useful for the calculation of the
    % stiffness block are stored.       
    LocEdge =  struct('id',id,'nini',Edges.nini(id), 'nfin',Edges.nfin(id),...
        'parametric',Edges.parametric(id,:),'lleft',Edges.lleft(id),...
        'lright',Edges.lright(id));
    
    %% Generating the geometric data
    % The following code transforms the abscissa coordinates, expressed in
    % the [-1,1] referential, to the polar coordinates required to compute
    % the values of the basis functions. The components of the outward
    % normal to the boundary in the radial and tangential directions are
    % also calculated. They are required to transform the stress basis into
    % the traction basis.
    
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

    % Computing the components of the outward normal in the polar
    % directions.
    NR = nx * cos(T) + ny * sin(T);   
    NT = -1*nx * sin(T) + ny * cos(T);
    
    
    %% Computing the basis functions for all integration points
    % Polar components of the displacement basis 
    Ur = R.^N.*sin((N+1).*T) ;
    Ut = R.^N.*cos((N+1).*T) ;
    % Polar components of the stress basis    
    Sr = 2*LocLoop.materials(5)*(LocLoop.materials(4)+...
        LocLoop.materials(5))*M.*R.^(M-1).*((M-3).*cos((M-1).*T));
    St = 2*LocLoop.materials(5)*(LocLoop.materials(4)+...
        LocLoop.materials(5))*M.*R.^(M-1).*(-(M+1).*cos((M-1).*T));
    Srt= 2*LocLoop.materials(5)*(LocLoop.materials(4)+...
        LocLoop.materials(5))*M.*R.^(M-1).*(-(M-1).*sin((M-1).*T));
    
    
    %% Computing the integral on the side
    % The integral is the internal product between the displacement basis U
    % and the traction basis N * S
    
    NUS = NR.*((Ur.*Sr)+(Ut.*Srt)) + NT.*((Ur.*Srt)+(Ut.*St));
    
    % Performing the side integration and updating the K14 matrix                                            
    w3D(1,1,:) = weight;
    K14i = K14i + L/2 * sum(bsxfun(@times,NUS,w3D),3); % computes the integral
 
end

end