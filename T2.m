function RHS = T2(Edges, Loops, BConds,RHS, abscissa, weight)
% T2 sweeps through the elements and calls the functions that generate the
% T2 blocks of the free vector (RHS) of the solving system. T2 blocks are
% only calculated for elements that contain Neumann boundaries. For all
% other elements, the T2 blocks are filled with zeros.
%
% T2 is called by MAIN***. It receives as input data the Edges, Loops and
% BConds structures, the RHS vector (that is, the free vector of the
% solving system), and the Gauss-Legendre integration parameters abscissa
% and weights. It returns to MAIN*** the RHS vector with the T2 blocks of
% all elements inserted at the correct positions (as determined in
% ASSIGNPARTS).
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
% T2 computes the internal product between the displacement basis U2 of the 
% element, expressed in a normal-tangential referential, and the applied 
% tractions on its Neumann boundaries. 
% * basis U2 is written in the polar (r,th) referential as,
%        U2(r,th) = | Ur | = r^n * | cos[(n+1)*th] |
%                   | Ut |         |-sin[(n+1)*th] | 
% where n is the line of the current entry in vector T2, and transformed to 
% the normal-tangential referential,
%        Un = nr*Ur + nt*Ut
%        Ut = -nt*Ur + nr*Ut
% where nr and nt are the radial and tangential components of the outward 
% unit normal to the current boundary.
% * the boundary traction function is defined in the GUI by its values in
% an arbitrary number of equally-spaced points along the boundary (see
% Section 4.6 of reference [3]). A polynomial interpolation is performed
% between these values to obtain the analytic expression of the applied
% tractions. The degree of the polynomial is equal to the number of
% traction values, minus one.
%
% Further details on the structure of the solving system are presented in 
% reference [2] (Section 6.2). The transformation of referentials is
% covered in reference [4] (Appendix C). 

%% Sweeping the elements
for ii=1:length(Loops.area)
    
    % LocLoop is a local structure where the features of the current
    % element which are directly useful for the calculation of the
    % T2 block are stored.
    LocLoop = struct('id',ii,'edges',Loops.edges(ii,:), 'center',...
        Loops.center(ii,:),'order',Loops.order(ii,1),... 
        'insert',Loops.insert(ii,2),'dim',Loops.dim(ii,2),...
        'materials',Loops.materials(ii,:));
    
    % Computing the T2 vector of element ii. Function T2_VECTOR_I is a
    % local function (see below).
    T2i = T2_Vector_i(LocLoop, Edges, BConds, abscissa, weight);
    
    % Inserting the T2 block in the RHS vector. The insertion is made at
    % line Loops.insert(ii,2).
    RHS(LocLoop.insert:LocLoop.insert+LocLoop.dim-1) = T2i;
    
end

end

function T2i = T2_Vector_i(LocLoop, Edges, BConds,abscissa, weight)
% T2_VECTOR_I local function computes the T2 vector of the LocLoop element. 
% The Neumann sides are mapped to a [-1,1] interval to perform the 
% integrations.

%% Initialization 
% Initialization of the T2 block
T2i = zeros(LocLoop.dim(1),1); 

% n denotes the current line
n = 1:LocLoop.order(1);

% Sweeping the edges for contour integration. Only Neumann boundaries will
% be selected to perform the internal product.
for jj = 1:length(LocLoop.edges)  
    
    % identification of the jj-th edge of the loop
    id = LocLoop.edges(jj); 
    
    % Neumann boundaries are selected to perform the internal product.
    if strcmpi(Edges.type(id),'N')    
        
        % LocEdge is a local structure where the features of the current
        % edge which are directly useful for the calculation of the
        % T2 block are stored.
        LocEdge =  struct('id',id,'nini',Edges.nini(id), 'nfin',Edges.nfin(id),...
            'parametric',Edges.parametric(id,:),'lleft',Edges.lleft(id),...
            'lright',Edges.lright(id));
        
        % Consistency check: exterior Neumann sides cannot have right loops
        if LocEdge.lright  
            error('local:consistencyChk',...
                'Exterior edge %d cannot have a right element. \n',...
                LocEdge.id);
        end
        
        %% Generating the geometric data
        % The following code transforms the abscissa coordinates, expressed 
        % in the [-1,1] referential, to the polar coordinates required to 
        % compute the values of the basis functions. The components of the 
        % outward normal to the boundary in the radial and tangential 
        % directions are also calculated. They are required to compute the 
        % normal and tangential components of the displacement basis.
        
        % Computing the length of the current edge
        L = sqrt(LocEdge.parametric(3)^2 + LocEdge.parametric(4)^2); 
        
        % Constructing the matrices containing the n x abscissa
        % integration grid
        [N,A] = ndgrid(n,abscissa);
        
        % Transforming the edge abscissa into local coordinates. The local
        % referential is centered in the barycenter of the element, its 
        % axes aligned with the Cartesian axes of the global referential.
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
        
        %% Computing the basis and traction functions at integration points
        % Polar components of the displacement basis
        Ur = R.^N.*cos((N+1).*T) ;
        Ut = R.^N.*-sin((N+1).*T) ;
        % Normal and tangential projections of the displacement basis
        Un = NR.*Ur + NT.*Ut;
        Utg = -NT.*Ur + NR.*Ut;
        
        % Computing the values of the applied tractions at the abscissas.
       
        % obtaining the equally spaced points on [-1,1] interval where the
        % tractions are defined and stored in BConds.Neumann
        an = linspace(-1,1,length(BConds.Neumann{id,1}));
        atg = linspace(-1,1,length(BConds.Neumann{id,2}));
        % obtaining the polynomials that interpolate the values in 
        % BConds.Neumann 
        pol_n = polyfit(an,BConds.Neumann{id,1},...
            length(BConds.Neumann{id,1})-1);
        pol_tg = polyfit(atg,BConds.Neumann{id,2},...
            length(BConds.Neumann{id,2})-1);
            
        % computing the values of the interpolation polynomials at the 
        % abscissas
        tn = polyval(pol_n,A);
        ttg = polyval(pol_tg,A);  
        
        %% Computing the integral on the side
        % The integral is the internal product between the displacement 
        % basis and the applied tractions in the normal and tangential
        % directions
        T2i2D = Un.*tn + Utg.*ttg; 
        
        % Performing the side integration and updating the T2 vector
        T2i = T2i + L/2 * sum(bsxfun(@times,T2i2D,weight.'),2);         
    end   
end
end