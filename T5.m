function RHS = T5(Edges, Loops, BConds,RHS, abscissa, weight)
% sweeps through the elements and calls the functions that generate the
% Xd vector in the RHS
% ************************************************************************
% Uses 1D and 2D data structures for storing n and the abscissas
% corresponding to all points/orders that must be computed. For integration
% it constructs a 3D matrix, with each page corresponding to the integrands 
% computed at a Gauss point. The integration is performed as a weighted 
% summation on the pages of this 3D matrix. 

for ii=1:length(Loops.area)
    
    LocLoop = struct('id',ii,'edges',Loops.edges(ii,:), 'center',...
        Loops.center(ii,:),'order',Loops.order(ii,1),... 
        'insert',Loops.insert(ii,5),'dim',Loops.dim(ii,5),...
        'materials',Loops.materials(ii,:));
    
    % Computing the Xd vector of element ii
    X1ri = X1r_Vector_i(LocLoop, Edges, BConds, abscissa, weight);
    
    % Inserting the vector in the global RHS vector
    RHS(LocLoop.insert) = X1ri;
    
        % Computing the Xd vector of element ii
    X2ri = X2r_Vector_i(LocLoop, Edges, BConds, abscissa, weight);
    
    % Inserting the vector in the global RHS vector
    RHS(LocLoop.insert+1) = X2ri;
    
        % Computing the Xd vector of element ii
    X3ri = X3r_Vector_i(LocLoop, Edges, BConds, abscissa, weight);
    
    % Inserting the vector in the global RHS vector
    RHS(LocLoop.insert+2) = X3ri;
    
end

end

function X1ri = X1r_Vector_i(LocLoop, Edges, BConds,abscissa, weight)

% computes the Xd vector of element ii

X1ri = 0;

n = 1;

for jj = 1:length(LocLoop.edges)  % contour integration
    
    id = LocLoop.edges(jj);  % number of the jj-th edge of the loop
    
    if strcmpi(Edges.type(id),'N')    
        
        LocEdge =  struct('id',id,'nini',Edges.nini(id), 'nfin',Edges.nfin(id),...
            'parametric',Edges.parametric(id,:),'lleft',Edges.lleft(id),...
            'lright',Edges.lright(id));
        
        if LocEdge.lright  % exterior Neumann sides cannot have right loops
            error('local:consistencyChk',...
                'Exterior edge %d cannot have a right element. \n',...
                LocEdge.id);
        end
        
        L = sqrt(LocEdge.parametric(3)^2 + LocEdge.parametric(4)^2); % length
        
        [N,A] = ndgrid(n,abscissa);
        
        % *****************************************************************
        % Getting the R, T, NR, NT for all Gauss points -- sendo que T é
        % teta, e tg para direcção tangente.
        
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
        Un = NR.*Ur + NT.*Ut;
        Utg = -NT.*Ur + NR.*Ut;
        
        
        % Computing the values of the flux at the abscissas:
       
        % obtaining the equally spaced points on [-1,1] interval where the
        % fluxes are defined and stored in BConds.Neumann
        
        an = linspace(-1,1,length(BConds.Neumann{id,1}));
        atg = linspace(-1,1,length(BConds.Neumann{id,2}));
        % obtaining the polynomial that gets the values in BConds.Neumann
        % at the points a
        
        pol_n = polyfit(an,BConds.Neumann{id,1},length(BConds.Neumann{id,1})-1);
        pol_tg = polyfit(atg,BConds.Neumann{id,2},length(BConds.Neumann{id,2})-1);
            
        % computing the values of "pol" at the abscissas
        
        tn = polyval(pol_n,A);
        ttg = polyval(pol_tg,A);        
%         tn = 1; ttg=1;
        X1ri2D = Un.*tn + Utg.*ttg; % this creates the 2D Xdi2D matrix,
                                            % one Gauss point per column
        
        X1ri = X1ri + L/2 * sum(bsxfun(@times,X1ri2D,weight.'),2); % computes the integral
        
    end
    
end

end


function X2ri = X2r_Vector_i(LocLoop, Edges, BConds,abscissa, weight)

% computes the Xd vector of element ii

X2ri = 0;

n = 1;

for jj = 1:length(LocLoop.edges)  % contour integration
    
    id = LocLoop.edges(jj);  % number of the jj-th edge of the loop
    
    if strcmpi(Edges.type(id),'N')    
        
        LocEdge =  struct('id',id,'nini',Edges.nini(id), 'nfin',Edges.nfin(id),...
            'parametric',Edges.parametric(id,:),'lleft',Edges.lleft(id),...
            'lright',Edges.lright(id));
        
        if LocEdge.lright  % exterior Neumann sides cannot have right loops
            error('local:consistencyChk',...
                'Exterior edge %d cannot have a right element. \n',...
                LocEdge.id);
        end
        
        L = sqrt(LocEdge.parametric(3)^2 + LocEdge.parametric(4)^2); % length
        
        [N,A] = ndgrid(n,abscissa);
        
        % *****************************************************************
        % Getting the R, T, NR, NT for all Gauss points -- sendo que T é
        % teta, e tg para direcção tangente.
        
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
                
        Ur = cos(T);
        Ut = -sin(T);
        Un = NR.*Ur + NT.*Ut;
        Utg = -NT.*Ur + NR.*Ut;
        
        
        % Computing the values of the flux at the abscissas:
       
        % obtaining the equally spaced points on [-1,1] interval where the
        % fluxes are defined and stored in BConds.Neumann
        
        an = linspace(-1,1,length(BConds.Neumann{id,1}));
        atg = linspace(-1,1,length(BConds.Neumann{id,2}));
        % obtaining the polynomial that gets the values in BConds.Neumann
        % at the points a
        
        pol_n = polyfit(an,BConds.Neumann{id,1},length(BConds.Neumann{id,1})-1);
        pol_tg = polyfit(atg,BConds.Neumann{id,2},length(BConds.Neumann{id,2})-1);
            
        % computing the values of "pol" at the abscissas
        
        tn = polyval(pol_n,A);
        ttg = polyval(pol_tg,A);        
%         tn = 1; ttg=1;
        X2ri2D = Un.*tn + Utg.*ttg; % this creates the 2D Xdi2D matrix,
                                            % one Gauss point per column
        
        X2ri = X2ri + L/2 * sum(bsxfun(@times,X2ri2D,weight.'),2); % computes the integral
        
    end
    
end

end


function X3ri = X3r_Vector_i(LocLoop, Edges, BConds,abscissa, weight)

% computes the Xd vector of element ii

X3ri = 0;

n = 1;

for jj = 1:length(LocLoop.edges)  % contour integration
    
    id = LocLoop.edges(jj);  % number of the jj-th edge of the loop
    
    if strcmpi(Edges.type(id),'N')    
        
        LocEdge =  struct('id',id,'nini',Edges.nini(id), 'nfin',Edges.nfin(id),...
            'parametric',Edges.parametric(id,:),'lleft',Edges.lleft(id),...
            'lright',Edges.lright(id));
        
        if LocEdge.lright  % exterior Neumann sides cannot have right loops
            error('local:consistencyChk',...
                'Exterior edge %d cannot have a right element. \n',...
                LocEdge.id);
        end
        
        L = sqrt(LocEdge.parametric(3)^2 + LocEdge.parametric(4)^2); % length
        
        [N,A] = ndgrid(n,abscissa);
        
        % *****************************************************************
        % Getting the R, T, NR, NT for all Gauss points -- sendo que T é
        % teta, e tg para direcção tangente.
        
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
                
        Ur = sin(T);
        Ut = cos(T);
        Un = NR.*Ur + NT.*Ut;
        Utg = -NT.*Ur + NR.*Ut;
        
        
        % Computing the values of the flux at the abscissas:
       
        % obtaining the equally spaced points on [-1,1] interval where the
        % fluxes are defined and stored in BConds.Neumann
        
        an = linspace(-1,1,length(BConds.Neumann{id,1}));
        atg = linspace(-1,1,length(BConds.Neumann{id,2}));
        % obtaining the polynomial that gets the values in BConds.Neumann
        % at the points a
        
        pol_n = polyfit(an,BConds.Neumann{id,1},length(BConds.Neumann{id,1})-1);
        pol_tg = polyfit(atg,BConds.Neumann{id,2},length(BConds.Neumann{id,2})-1);
            
        % computing the values of "pol" at the abscissas
        
        tn = polyval(pol_n,A);
        ttg = polyval(pol_tg,A);        
%         tn = 1; ttg=1;
        X3ri2D = Un.*tn + Utg.*ttg; % this creates the 2D Xdi2D matrix,
                                            % one Gauss point per column
        
        X3ri = X3ri + L/2 * sum(bsxfun(@times,X3ri2D,weight.'),2); % computes the integral
        
    end
    
end

end