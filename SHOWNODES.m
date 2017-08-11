function SHOWNODES(ShowNodes,empty,ShowEdges,ShowElements,nodes,edges,loops_nodes,X,Y,nnode,nedge,nel) 
%--------------------------------------------------------------------------
%   Purpose:
%           To display and undisplay the node numbers
%--------------------------------------------------------------------------
ShowElements = get(ShowElements,'Value') ;
ShowEdges = get(ShowEdges,'Value') ;
ShowNodes = get(ShowNodes,'Value') ;  

delta = sqrt(max(max(X))^2+max(max(Y))^2);  % to scale the insertion points of the text

% Display only Element Edges

cla(gcf) ;
patch(X,Y,'w') ;

if ShowNodes==1
    for i = 1:nnode
        text(nodes(i,1)+0.0*delta,nodes(i,2)+0.0*delta,int2str(i),'fontsize',8,....
            'fontweight','bold','Color','r'); 
    end    
end

if ShowEdges==1   
    for i = 1:nedge 
        EX = nodes(edges(i,:),1) ; EY = nodes(edges(i,:),2) ;
        pos = [sum(EX)/2+0.0*delta,sum(EY)/2+0.0*delta] ;
        text(pos(1),pos(2),int2str(i),'fontsize',8, ...
            'fontweight','bold','color','g');
    end
end

if ShowElements==1 
    for i = 1:nel
        C = polygonCentroid(nodes(loops_nodes(i,:),:));
        pos = C ;
        text(pos(1),pos(2),int2str(i),'fontsize',8, 'fontweight','bold',...
            'color','b');
    end
end

return