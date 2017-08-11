function RHS = U(Edges, BConds, RHS, abscissa, weight)
% sweeps through the edges and calls the functions that generate the
% Q vectors on the exterior Dirichlet sides
% ************************************************************************
% Uses 1D data structures for storing n and abscissas
% corresponding to all points/orders that must be computed. For integration
% it constructs a 2D matrix, with each column corresponding to the integrand 
% computed at a Gauss point. The integration is performed as a weighted 
% summation on the columns of this 2D matrix. 

for ii=1:length(Edges.type)
    
    if (strcmpi(Edges.type(ii),'D') && Edges.lright(ii) == 0 )
        LocEdge =struct('id',ii,'nini',Edges.nini(ii),'nfin',Edges.nfin(ii),...
            'parametric',Edges.parametric(ii,:), 'lleft',Edges.lleft(ii),...
            'lright',Edges.lright(ii), 'order',Edges.order(ii),...
            'insert',Edges.insert(ii,:), 'dim',Edges.dim(ii,:));
        
        for jj = 1:2 % normal & tangential
            
            if LocEdge.dim(jj)
                
                % Computing the Q vector of edge ii
                Qi = Q_Vector_i(LocEdge, BConds, jj, abscissa, weight);
                
                % Inserting the vector in the global RHS vector
                RHS(LocEdge.insert(jj):LocEdge.insert(jj)+LocEdge.dim(jj)-1) = -Qi;
            end
        end
        
    end
    
end
end

function Qi = Q_Vector_i(LocEdge, BConds, jj, abscissa, weight)

% computes the Q vector for edge LocEdge 

Qi = zeros(LocEdge.dim(jj),1);

Qi2D = zeros(LocEdge.dim(jj),length(abscissa));

m = 0:LocEdge.order;

L = sqrt(LocEdge.parametric(3)^2 + LocEdge.parametric(4)^2); % length

% Integrating on the side (the side integration is fully vectorialized)

% Z* -> the order is 'n'
Z = cos(bsxfun(@times,m,acos(abscissa)));
% Z = Z.';

% Computing the values of the enforced "displacements" at the abscissas:

% obtaining the equally spaced points on [-1,1] interval where the
% "displacements" are defined and stored in BConds.Dirichlet

% a = linspace(-1,1,length(BConds.Dirichlet{LocEdge.id}));
a = linspace(-1,1,length(BConds.Dirichlet{LocEdge.id,jj}));


% obtaining the polynomial that has the values given in BConds.Dirichlet
% at the points a

if (isnan(BConds.Dirichlet{LocEdge.id,jj}))
    error('local:consistencyChk',...
        'No Dirichlet boundary conditions are defined on edge %d. \n',...
        LocEdge.id(jj));
else
    pol = polyfit(a,BConds.Dirichlet{LocEdge.id,jj},...
        length(BConds.Dirichlet{LocEdge.id,jj})-1);
end

% computing the values of "pol" at the abscissas

q = polyval(pol,abscissa);
% q = q.'; % non-conjugate transpose
% q = 18;
Qi2D = bsxfun(@times, Z, q); % this creates the 2D Qi2D matrix,
% one Gauss point per column

Qi = L/2 * sum(bsxfun(@times,Qi2D,weight),1); % computes the integral
Qi = Qi.';

end