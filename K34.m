function LHS = K34(Edges, Loops, LHS, abscissa, weight)
% sweeps through the elements and calls the functions that generate the
% K34 stiffness matrix in the LHS.

for ii=1:length(Loops.area)
    
    LocLoop = struct('id',ii,'edges',Loops.edges(ii,:), 'center',...
        Loops.center(ii,:),'order',Loops.order(ii,:),... 
        'insert',Loops.insert(ii,:),'dim',Loops.dim(ii,:),...
        'materials',Loops.materials(ii,:));
    
    % Computing the Ddd matrix of element ii
    K34i = K34_Matrix_i(LocLoop, Edges, abscissa, weight);
    
    % Inserting the matrix in the global LHS matrix
    LHS(LocLoop.insert(3):LocLoop.insert(3)+LocLoop.dim(3)-1,...
        LocLoop.insert(4):LocLoop.insert(4)+LocLoop.dim(4)-1) = K34i;
    % Inserting the transposed of the matrix in the global LHS matrix
    LHS(LocLoop.insert(4):LocLoop.insert(4)+LocLoop.dim(4)-1,...
        LocLoop.insert(3):LocLoop.insert(3)+LocLoop.dim(3)-1) = K34i';
    
end

end

function K34i = K34_Matrix_i(LocLoop, Edges, abscissa, weight)

% computes the Ddd matrix of element ii

K34i = zeros([LocLoop.dim(3),LocLoop.dim(4)]);

n = 2:LocLoop.order(1);
m = 1:LocLoop.order(1);

for jj = 1:length(LocLoop.edges)  % contour integration
    
    id = LocLoop.edges(jj);  % number of the jj-th edge of the loop
    LocEdge =  struct('id',id,'nini',Edges.nini(id), 'nfin',Edges.nfin(id),...
        'parametric',Edges.parametric(id,:),'lleft',Edges.lleft(id),...
        'lright',Edges.lright(id));
    
    L = sqrt(LocEdge.parametric(3)^2 + LocEdge.parametric(4)^2); % length
    
    % Constructing the 3D matrices containing the n x m x abscissa
    % integration grid
    
    [N,M,A] = ndgrid(n,m,abscissa);
    
    % *****************************************************************
    % Getting the r, th, nr, nth for all Gauss points
    
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
    
    NR = nx * cos(T) + ny * sin(T);   % normal in local r-th
    NT = -1*nx * sin(T) + ny * cos(T);
    
    % *****************************************************************

    
        
   % base Uhs e Shs (44);
    Ur = R.^N.*sin((N-1).*T).*((LocLoop.materials(5)*(N-3))+(LocLoop.materials(4)*(N-1)));
    Ut = R.^N.*cos((N-1).*T).*((LocLoop.materials(5)*(N+3))+(LocLoop.materials(4)*(N+1))) ;
    Sr = 2*LocLoop.materials(5)*(LocLoop.materials(4)+LocLoop.materials(5))*M.*R.^(M-1).*((M-3).*cos((M-1).*T));
    St = 2*LocLoop.materials(5)*(LocLoop.materials(4)+LocLoop.materials(5))*M.*R.^(M-1).*(-(M+1).*cos((M-1).*T));
    Srt= 2*LocLoop.materials(5)*(LocLoop.materials(4)+LocLoop.materials(5))*M.*R.^(M-1).*(-(M-1).*sin((M-1).*T));
    
    
    %Base de aproximção, que se integrará;
    
    NUS = NR.*((Ur.*Sr)+(Ut.*Srt)) + NT.*((Ur.*Srt)+(Ut.*St));
                                            
    w3D(1,1,:) = weight;
    
    K34i = K34i + L/2 * sum(bsxfun(@times,NUS,w3D),3); % computes the integral
 
end

end