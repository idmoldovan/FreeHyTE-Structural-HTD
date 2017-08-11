function [nodes, edges_nodes, edges_loops, loops_nodes, loops_edges] = ... 
    CreateEdgeLoop(p,t)

nodes = p';

t(end,:)=t(1,:); % drops the subdomain line of t and substitutes it with the 
               % first, to facilitate the construction of edge structures

% dumpedge contains the nodes of all edges, even redundant. The third line
% contains the respective element.
dumpedge = zeros(3,(size(t,1)-1)*size(t,2)); % '-1' as the last line is redundant
k = 1;
for i = 1:size(t,2) % number of loops
    for j = 1:size(t,1)-1 % number of nodes per loop is equal to the number of edges
        dumpedge(1:2,k) = t(j:j+1,i);
        dumpedge(3,k) = i;
        k=k+1;
    end
end
dumpedge = dumpedge';

% start the filtering process on dumpedge
edges_nodes = [];
edges_loops = [];
k=1;
for i = 1:size(dumpedge,1)
   crtedge = dumpedge(i,1:2); 
   
   % checking if crtedge was already registered in edge_nodes
   if  find(sum(ismember(edges_nodes,crtedge),2)>1)
       % sum(ismember(edges_nodes,crtedge),2) checks for the elements of crtedge
       % in edges_nodes and yields a vector containing 0 if none of the
       % elements of crtedge is present on the respective line of
       % edges_nodes, 1 if one is present and 2 if both are present. Then,
       % we are looking to see if there are lines with value 2 - this means
       % that crtedge is present and we want to skip it.
       
       continue;
   else % if the edge is not present, we register it
       edges_nodes(k,:) = crtedge;
       
       edges_loops(k,1) = dumpedge(i,3); % the left element
       
       % identifying other locations in dumpedge where crtedge occurs.
       % setdiff looks for other locations than 'i'.
       diffloc = setdiff(find(sum(ismember(dumpedge(:,1:2),crtedge),2)>1),i);
       if diffloc % if there are other locations where crtedge occurs
           edges_loops(k,2) = dumpedge(diffloc,3); % the right element
       else
           edges_loops(k,2) = 0; % no right element
       end
       k = k+1;
   end
    
end

% construct loops_edges: take dumpedge and for each entry identify the
% corresponding edge in edges_nodes. Take its number and store it into the
% loops_edges matrix.
loops_edges = zeros(size(t,2),size(t,1)-1); % contains as many lines as elements and as many columns as the elements' nodes.

for i = 1:size(dumpedge,1)
    crtedge = dumpedge(i,1:2);
    edgenumber = find(sum(ismember(edges_nodes,crtedge),2)>1);
    loopnumber = dumpedge(i,3);
    
    % edgenumber must now be inserted into the loops_edges structure, at
    % line loopnumber and in the first column with an empty (zero) entry.
    
    % finding the next zero entry on lineloopnumber
    colnumber = find(loops_edges(loopnumber,:)==0, 1 );
    
    loops_edges(loopnumber,colnumber) = edgenumber;
    
end

loops_nodes = t(1:3,:)';

end