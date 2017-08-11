function LHS = B5(Edges, Loops, LHS, abscissa, weight)
% sweeps through the edges and calls the functions that generate the
% Bs matrices for the left and right loops
% ************************************************************************
% Uses 1D and 2D data structures for storing n, m and abscissas
% corresponding to all points/orders that must be computed. For integration
% it constructs a 3D matrix, with each page corresponding to the integrands 
% computed at a Gauss point. The integration is performed as a weighted 
% summation on the pages of this 3D matrix. 

for ii=1:length(Edges.type)
    
    if strcmpi(Edges.type(ii),'D')
        LocEdge =struct('nini',Edges.nini(ii),'nfin',Edges.nfin(ii),...
            'parametric',Edges.parametric(ii,:), 'lleft',Edges.lleft(ii),...
            'lright',Edges.lright(ii), 'order',Edges.order(ii),...
            'insert',Edges.insert(ii,:), 'dim',Edges.dim(ii,:));
        
        if LocEdge.lleft
            id = LocEdge.lleft;
            sign = 1.;
            LocLoop = struct('id',id,'edges',Loops.edges(id,:),... 
                'center',Loops.center(id,:),'order',Loops.order(id,1),...
                'insert',Loops.insert(id,5),'dim',Loops.dim(id,5),...
                'materials',Loops.materials(id,:));
            
            for jj = 1:2 % normal & tangential
                
                if LocEdge.dim(jj)
                   
                    % Computing the Bs matrix of edge ii, left loop
                    B1ri = sign*B1r_Matrix_i(LocEdge, LocLoop, jj, abscissa, weight);

                    % Inserting the matrix in the global LHS matrix
                    LHS(LocLoop.insert,...
                        LocEdge.insert(jj):LocEdge.insert(jj)+LocEdge.dim(jj)-1) = -B1ri;
                    % Inserting the conjugate transposed in the global LHS matrix
                    LHS(LocEdge.insert(jj):LocEdge.insert(jj)+LocEdge.dim(jj)-1,...
                        LocLoop.insert) = -B1ri';
                    
                   
                    
                    
                    % Computing the Bs matrix of edge ii, left loop
                    B2ri = sign*B2r_Matrix_i(LocEdge, LocLoop, jj, abscissa, weight);
                    
                    % Inserting the matrix in the global LHS matrix
                    LHS(LocLoop.insert +1 ,...
                        LocEdge.insert(jj):LocEdge.insert(jj)+LocEdge.dim(jj)-1) = -B2ri;
                    % Inserting the conjugate transposed in the global LHS matrix
                    LHS(LocEdge.insert(jj):LocEdge.insert(jj)+LocEdge.dim(jj)-1,...
                        LocLoop.insert+1 ) = -B2ri';
                    
                    
                    
                    
                    % Computing the Bs matrix of edge ii, left loop
                    B3ri = sign*B3r_Matrix_i(LocEdge, LocLoop, jj, abscissa, weight);
                    
                    % Inserting the matrix in the global LHS matrix
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
        
        if LocEdge.lright
            id = LocEdge.lright;
            sign = +1.;
            LocLoop = struct('id',id,'edges',Loops.edges(id,:),... 
                'center',Loops.center(id,:),'order',Loops.order(id,1),...
                'insert',Loops.insert(id,5),'dim',Loops.dim(id,5));
            
            for jj = 1:2 % normal & tangential
                
                if LocEdge.dim(jj)
               
                    % Computing the Bs matrix of edge ii, right loop
                    B1ri = sign*B1r_Matrix_i(LocEdge, LocLoop, jj, abscissa, weight);
            
                    % Inserting the matrix in the global LHS matrix
                    LHS(LocLoop.insert,...
                        LocEdge.insert(jj):LocEdge.insert(jj)+LocEdge.dim(jj)-1) = -B1ri;
                    % Inserting the conjugate transposed in the global LHS matrix
                    LHS(LocEdge.insert(jj):LocEdge.insert(jj)+LocEdge.dim(jj)-1,...
                        LocLoop.insert) = -B1ri';
                    
                    
                    
                    
                    % Computing the Bs matrix of edge ii, right loop
                    B2ri = sign*B2r_Matrix_i(LocEdge, LocLoop, jj, abscissa, weight);
                    
                    % Inserting the matrix in the global LHS matrix
                    LHS(LocLoop.insert+1,...
                        LocEdge.insert(jj):LocEdge.insert(jj)+LocEdge.dim(jj)-1) = -B2ri;
                    % Inserting the conjugate transposed in the global LHS matrix
                    LHS(LocEdge.insert(jj):LocEdge.insert(jj)+LocEdge.dim(jj)-1,...
                        LocLoop.insert+1) = -B2ri';
                    
                    
                    
                    
                    % Computing the Bs matrix of edge ii, right loop
                    B3ri = sign*B3r_Matrix_i(LocEdge, LocLoop, jj, abscissa, weight);
                    
                    % Inserting the matrix in the global LHS matrix
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

% computes the Bs matrix for edge LocEdge and loop LocLoop

B1ri = zeros(LocLoop.dim,LocEdge.dim(jj));

B1ri3D = zeros(LocLoop.dim,LocEdge.dim(jj),length(abscissa));

% Z = zeros(LocEdge.dim,LocEdge.dim);

% n = 1:LocLoop.order(1);
n = 1;
m = 0:LocEdge.order(1);                      

L = sqrt(LocEdge.parametric(3)^2 + LocEdge.parametric(4)^2); % length

[N,M,A] = ndgrid(n,m,abscissa);



% *****************************************************************
% Getting the R, T, NR, NT for all Gauss points

loc_x = LocEdge.parametric(1) - LocLoop.center(1) + 0.5 *...
    (A + 1) * LocEdge.parametric(3);  % x & y in local ccord
loc_y = LocEdge.parametric(2) - LocLoop.center(2) + 0.5 *...
    (A + 1) * LocEdge.parametric(4);

R = sqrt(loc_x.^2 + loc_y.^2);  % polar coordinates, local
T = atan2(loc_y, loc_x);

nx = LocEdge.parametric(4) / L;   % normal in (local/global) x & y
        ny = -1* LocEdge.parametric(3) / L;
        if LocEdge.lright==LocLoop.id  % if the element is on the right,
            nx = -nx;                 % change the sign of the normal
            ny = -ny;
        end
        
        NR = nx * cos(T) + ny * sin(T);   % normal in local R-T
        NT = -1*nx * sin(T) + ny * cos(T);
 
        
        
% *****************************************************************
% Integrating on the side (the side integration is fully vectorialized)

% U* -> the order is 'n'
Ur = 0;
Ut = R;
if jj == 1
    U = NR.*Ur + NT.*Ut;
else
    U = -NT.*Ur + NR.*Ut;
end

% z -> Chebyshev functions of degree 'm'
Z = cos(bsxfun(@times,M,acos(A)));

% Z = [z 0;0 z];

B1ri3D = U.*Z;
% B1i3D = bsxfun(@times, Ustar3D, Z3D); % this creates the 3D Bsi matrix,
% one Gauss point per page

w3D(1,1,:) = weight;

B1ri = L/2 * sum(bsxfun(@times,B1ri3D,w3D),3); % computes the integral


end



function B2ri = B2r_Matrix_i(LocEdge, LocLoop, jj, abscissa, weight)

% computes the Bs matrix for edge LocEdge and loop LocLoop

B2ri = zeros(LocLoop.dim,LocEdge.dim(jj));

B2ri3D = zeros(LocLoop.dim,LocEdge.dim(jj),length(abscissa));

% Z = zeros(LocEdge.dim,LocEdge.dim);

n = 1;
m = 0:LocEdge.order(1);                      

L = sqrt(LocEdge.parametric(3)^2 + LocEdge.parametric(4)^2); % length

[N,M,A] = ndgrid(n,m,abscissa);



% *****************************************************************
% Getting the R, T, NR, NT for all Gauss points

loc_x = LocEdge.parametric(1) - LocLoop.center(1) + 0.5 *...
    (A + 1) * LocEdge.parametric(3);  % x & y in local ccord
loc_y = LocEdge.parametric(2) - LocLoop.center(2) + 0.5 *...
    (A + 1) * LocEdge.parametric(4);

% R = sqrt(loc_x.^2 + loc_y.^2);  % polar coordinates, local
T = atan2(loc_y, loc_x);

nx = LocEdge.parametric(4) / L;   % normal in (local/global) x & y
        ny = -1* LocEdge.parametric(3) / L;
        if LocEdge.lright==LocLoop.id  % if the element is on the right,
            nx = -nx;                 % change the sign of the normal
            ny = -ny;
        end
        
        NR = nx * cos(T) + ny * sin(T);   % normal in local R-T
        NT = -1*nx * sin(T) + ny * cos(T);

        
        
% *****************************************************************
% Integrating on the side (the side integration is fully vectorialized)

% U* -> the order is 'n'
Ur = cos(T);
Ut = -sin(T);
if jj == 1
    U = NR.*Ur + NT.*Ut;
else
    U = -NT.*Ur + NR.*Ut;
end

% z -> Chebyshev functions of degree 'm'
Z = cos(bsxfun(@times,M,acos(A)));

% Z = [z 0;0 z];

B2ri3D = U.*Z;
% B1i3D = bsxfun(@times, Ustar3D, Z3D); % this creates the 3D Bsi matrix,
% one Gauss point per page

w3D(1,1,:) = weight;

B2ri = L/2 * sum(bsxfun(@times,B2ri3D,w3D),3); % computes the integral


end




function B3ri = B3r_Matrix_i(LocEdge, LocLoop, jj, abscissa, weight)

% computes the Bs matrix for edge LocEdge and loop LocLoop

B3ri = zeros(LocLoop.dim,LocEdge.dim(jj));

B3ri3D = zeros(LocLoop.dim,LocEdge.dim(jj),length(abscissa));

% Z = zeros(LocEdge.dim,LocEdge.dim);

n = 1;
m = 0:LocEdge.order(1);                      

L = sqrt(LocEdge.parametric(3)^2 + LocEdge.parametric(4)^2); % length

[N,M,A] = ndgrid(n,m,abscissa);



% *****************************************************************
% Getting the R, T, NR, NT for all Gauss points

loc_x = LocEdge.parametric(1) - LocLoop.center(1) + 0.5 *...
    (A + 1) * LocEdge.parametric(3);  % x & y in local ccord
loc_y = LocEdge.parametric(2) - LocLoop.center(2) + 0.5 *...
    (A + 1) * LocEdge.parametric(4);

% R = sqrt(loc_x.^2 + loc_y.^2);  % polar coordinates, local
T = atan2(loc_y, loc_x);

nx = LocEdge.parametric(4) / L;   % normal in (local/global) x & y
        ny = -1* LocEdge.parametric(3) / L;
        if LocEdge.lright==LocLoop.id  % if the element is on the right,
            nx = -nx;                 % change the sign of the normal
            ny = -ny;
        end
        
        NR = nx * cos(T) + ny * sin(T);   % normal in local R-T
        NT = -1*nx * sin(T) + ny * cos(T);

        
        
% *****************************************************************
% Integrating on the side (the side integration is fully vectorialized)

% U* -> the order is 'n'
Ur = sin(T);
Ut = cos(T);
if jj == 1
    U = NR.*Ur + NT.*Ut;
else
    U = -NT.*Ur + NR.*Ut;
end

% z -> Chebyshev functions of degree 'm'
Z = cos(bsxfun(@times,M,acos(A)));

% Z = [z 0;0 z];

B3ri3D = U.*Z;
% B1i3D = bsxfun(@times, Ustar3D, Z3D); % this creates the 3D Bsi matrix,
% one Gauss point per page

w3D(1,1,:) = weight;

B3ri = L/2 * sum(bsxfun(@times,B3ri3D,w3D),3); % computes the integral


end