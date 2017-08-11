function [Edges,Loops,Dim] = AssignParts(Edges, Loops, BConds)

%The point of this function is to assign to each element and Dirichlet
%side an entry point and a dimension in the global solving system
%--------------------------------------------------------------------------

entry = 1;

% initialize the insertion points for the static and dynamic parts
Loops.insert = zeros(length(Loops.area),5); 
% initialize the dimensions of the static and dynamic parts
Loops.dim = zeros(length(Loops.area),5);

for i = 1:length(Loops.area)
    Loops.insert(i,1) = entry;
    Loops.dim(i,1) = Loops.order(i);
    entry = entry + Loops.dim(i,1);
    Loops.insert(i,2) = entry;
    Loops.dim(i,2) = Loops.order(i);
    entry = entry + Loops.dim(i,2);
    Loops.insert(i,3) = entry;
    Loops.dim(i,3) = Loops.order(i)-1;
    entry = entry + Loops.dim(i,3);
    Loops.insert(i,4) = entry;
    Loops.dim(i,4) = Loops.order(i);
    entry = entry + Loops.dim(i,4);
    Loops.insert(i,5) = entry;
    Loops.dim(i,5) = 3;
    entry = entry + Loops.dim(i,5);
    
    
    
end

% initialize the insertion points for the edges
Edges.insert = zeros(length(Edges.type),2); 
% initialize the dimensions of the static and dynamic parts
Edges.dim = zeros(length(Edges.type),2); 

for i = 1:length(Edges.insert(:,1))
    if strcmpi(Edges.type(i),'D')
        if(~Edges.lright(i))    % fronteira exterior
            if any(~isnan(BConds.Dirichlet{i,1})) % if any of the terms in BConds.Dirichlet{i,1} is different from NaN, it executes
                Edges.insert(i,1) = entry;
                Edges.dim(i,1) = Edges.order(i)+1;
                entry = entry + Edges.dim(i,1);
            end
            if any(~isnan(BConds.Dirichlet{i,2}))
                Edges.insert(i,2) = entry;
                Edges.dim(i,2) = Edges.order(i)+1;
                entry = entry + Edges.dim(i,2);
            end
        else     % fronteira interior
                Edges.insert(i,1) = entry;
                Edges.dim(i,1) = Edges.order(i)+1;
                entry = entry + Edges.dim(i,1);
                Edges.insert(i,2) = entry;
                Edges.dim(i,2) = Edges.order(i)+1;
                entry = entry + Edges.dim(i,2);
        end
    end
end  % note that Neumann edges correspond to zero dim, zero insert point

Dim = entry-1;

end

