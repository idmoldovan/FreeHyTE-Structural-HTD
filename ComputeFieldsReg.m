function [Ux,Uy,Sx,Sy,Sxy] = ComputeFieldsReg(NoDiv, Nodes, Edges, Loops, X)
% Computes and plots the final displacement and stress fields

% Creating the sub-folder to store the analysis results into
load('StructDef','DirName','FileName','EdgesOrder','LoopsOrder'); %loads the folder defined in the GUI

if ~isempty(DirName)  % if the user left DirName blank, does not generate written output
    
    % Generating the file name
    % removing .mat extension
    FileName = FileName(1:end-4);
    FileName = sprintf('%s_ND%dNB%d.dat',FileName,LoopsOrder,EdgesOrder);
    % link the path to the filename
    UFilename = fullfile(DirName,FileName);
    % open file for writing
    FileU = fopen(UFilename,'w');
    fprintf(FileU,'TITLE="%s"\n',FileName);
    fprintf(FileU,'VARIABLES="X", "Y", "Ux", "Uy", "Sx", "Sy", "Sxy" \n'); % \n functions as <ENT>
    fprintf(FileU,'ZONE T="ND = %d; NB = %d"\n',LoopsOrder,EdgesOrder);
    
end

% get the coordinates of the mesh points
xmesh = [Nodes(Edges.nini(:),1) Nodes(Edges.nfin(:),1)];
ymesh = [Nodes(Edges.nini(:),2) Nodes(Edges.nfin(:),2)];

% draw the mesh, once for Ux, once for Uy
Fig=figure;
set(Fig,'name','Displacement and Stress field','numbertitle','off','color','w') ;

Ux_Plt = subplot(3,2,1);
hold on; title('X- Displacement');
axis([min(min(xmesh)) max(max(xmesh)) min(min(ymesh)) max(max(ymesh))]);
for ii=1:length(Edges.type)
    line(xmesh(ii,:),ymesh(ii,:),'color','b');
end
daspect([1 1 1]);

Uy_Plt = subplot(3,2,2);
hold on; title('Y- Displacement');
axis([min(min(xmesh)) max(max(xmesh)) min(min(ymesh)) max(max(ymesh))]);
for ii=1:length(Edges.type)
    line(xmesh(ii,:),ymesh(ii,:),'color','b');
end
daspect([1 1 1]);

Sx_Plt = subplot(3,2,3);
hold on; title('X- Stress');
axis([min(min(xmesh)) max(max(xmesh)) min(min(ymesh)) max(max(ymesh))]);
for ii=1:length(Edges.type)
    line(xmesh(ii,:),ymesh(ii,:),'color','b');
end
daspect([1 1 1]);

Sy_Plt = subplot(3,2,4);
hold on; title('Y- Stress');
axis([min(min(xmesh)) max(max(xmesh)) min(min(ymesh)) max(max(ymesh))]);
for ii=1:length(Edges.type)
    line(xmesh(ii,:),ymesh(ii,:),'color','b');
end
daspect([1 1 1]);

Sxy_Plt = subplot(3,2,5);
hold on; title('XY- Stress');
axis([min(min(xmesh)) max(max(xmesh)) min(min(ymesh)) max(max(ymesh))]);
for ii=1:length(Edges.type)
    line(xmesh(ii,:),ymesh(ii,:),'color','b');
end
daspect([1 1 1]);

Def_Plt = subplot(3,2,6);
hold on; title('Deformed shape');
%axis([min(min(xmesh)) max(max(xmesh)) min(min(ymesh)) max(max(ymesh))]);
for ii=1:length(Edges.type)
    line(xmesh(ii,:),ymesh(ii,:),'color',[0.5,0.5,0.5]);
end
daspect([1 1 1]);

% ********************************************************************

% Generating the mesh of output points. A very small shift is induced to
% the limites of the elements in order to avoid the duplication of boundary
% points on adjacent elements, which may compromise the visibility of the
% continuity gaps.
pts = linspace(-(1-1.e-5),(1-1.e-5),NoDiv+1);

% Store the maximum displacement (in either direction) in the structure,
% to draw its deformed shape
Umax = 0;

% Computing & drawing the displacement and stress fields

for ii=1:length(Loops.area)
    
    LocLoop = struct('id',ii,'nodes',Loops.nodes(ii,:),'center',...
        Loops.center(ii,:),'order',Loops.order(ii),'insert',...
        Loops.insert(ii,:),'dim',Loops.dim(ii,:),...
        'materials',Loops.materials(ii,:));
    
    % Static modes x and y(3D, on pages)
    n = 1:LocLoop.order;
    
    
    % Computing the length of the sides of the element in x and y direction.
    % size simply collects the distances between the two most far apart points
    % of the element in x and y directions. ASSUMES THAT THE ELEMENT IS
    % RECTANGULAR!!
    sze = max(Nodes(LocLoop.nodes(:),:))-min(Nodes(LocLoop.nodes(:),:));
    
    % Generating the output grid (local coordinates)
    [x,y,N] = ndgrid(1/2*sze(1)*pts,1/2*sze(2)*pts,n);
    
    % Calculating the R, T coordinates
    R = sqrt(x.^2 + y.^2);
    T = atan2(y,x);
    
    % Modos de deslocamentos e tensões
    U1r = R.^N.*sin((N+1).*T);
    U1t = R.^N.*cos((N+1).*T);
    U2r = R.^N.*cos((N+1).*T);
    U2t = R.^N.*-sin((N+1).*T);
    U3r = R.^N.*sin((N-1).*T).*((LocLoop.materials(5)*(N-3))+(LocLoop.materials(4)*(N-1)));
    U3t = R.^N.*cos((N-1).*T).*((LocLoop.materials(5)*(N+3))+(LocLoop.materials(4)*(N+1)));
    U4r = R.^N.*cos((N-1).*T).*((LocLoop.materials(5)*(N-3))+(LocLoop.materials(4)*(N-1)));
    U4t = -R.^N.*sin((N-1).*T).*((LocLoop.materials(5)*(N+3))+(LocLoop.materials(4)*(N+1)));
    U1rrb = 0;
    U1trb = R;
    U2rrb = cos(T);
    U2trb = -sin(T);
    U3rrb = sin(T);
    U3trb = cos(T);
    
    S1r =2*LocLoop.materials(5)*N.*R.^(N-1).*sin((N+1).*T);
    S1t = 2*LocLoop.materials(5)*N.*R.^(N-1).*-sin((N+1).*T);
    S1rt = 2*LocLoop.materials(5)*N.*R.^(N-1).*cos((N+1).*T);
    S2r = 2*LocLoop.materials(5)*N.*R.^(N-1).*cos((N+1).*T);
    S2t = 2*LocLoop.materials(5)*N.*R.^(N-1).*-cos((N+1).*T);
    S2rt = 2*LocLoop.materials(5)*N.*R.^(N-1).*-sin((N+1).*T);
    S3r = 2*LocLoop.materials(5)*(LocLoop.materials(4)+LocLoop.materials(5))*N.*R.^(N-1).*((N-3).*sin((N-1).*T));
    S3t = 2*LocLoop.materials(5)*(LocLoop.materials(4)+LocLoop.materials(5))*N.*R.^(N-1).*(-(N+1).*sin((N-1).*T));
    S3rt = 2*LocLoop.materials(5)*(LocLoop.materials(4)+LocLoop.materials(5))*N.*R.^(N-1).*((N-1).*cos((N-1).*T));
    S4r = 2*LocLoop.materials(5)*(LocLoop.materials(4)+LocLoop.materials(5))*N.*R.^(N-1).*((N-3).*cos((N-1).*T));
    S4t = 2*LocLoop.materials(5)*(LocLoop.materials(4)+LocLoop.materials(5))*N.*R.^(N-1).*(-(N+1).*cos((N-1).*T));
    S4rt= 2*LocLoop.materials(5)*(LocLoop.materials(4)+LocLoop.materials(5))*N.*R.^(N-1).*(-(N-1).*sin((N-1).*T));
    
    % Vector X que multiplica com os modos (U e S)
    X1(1,1,:) = X(LocLoop.insert(1):LocLoop.insert(1)+LocLoop.dim(1)-1);
    X2(1,1,:) = X(LocLoop.insert(2):LocLoop.insert(2)+LocLoop.dim(2)-1);
    X3(1,1,:) = X(LocLoop.insert(3):LocLoop.insert(3)+LocLoop.dim(3)-1);
    X4(1,1,:) = X(LocLoop.insert(4):LocLoop.insert(4)+LocLoop.dim(4)-1);
    X5 = X(LocLoop.insert(5):LocLoop.insert(5)+LocLoop.dim(5)-1);
    
    % The first page in U3 is cut off from the basis (corresponds to a zero mode)
    Ur = sum(bsxfun(@times,U1r,X1),3) + sum(bsxfun(@times,U2r,X2),3) + ...
        sum(bsxfun(@times,U3r(:,:,2:end),X3),3) + sum(bsxfun(@times,U4r,X4),3)...
        + U1rrb(:,:,1)*X5(1) + U2rrb(:,:,1)*X5(2) + U3rrb(:,:,1)*X5(3); %
    
    Ut = sum(bsxfun(@times,U1t,X1),3) + sum(bsxfun(@times,U2t,X2),3) + ...
        sum(bsxfun(@times,U3t(:,:,2:end),X3),3) + sum(bsxfun(@times,U4t,X4),3)...
        + U1trb(:,:,1)*X5(1) + U2trb(:,:,1)*X5(2) + U3trb(:,:,1)*X5(3); %
    
    Sr = sum(bsxfun(@times,S1r,X1),3) + sum(bsxfun(@times,S2r,X2),3) + ...
        sum(bsxfun(@times,S3r(:,:,2:end),X3),3) + sum(bsxfun(@times,S4r,X4),3); %
    
    St = sum(bsxfun(@times,S1t,X1),3) + sum(bsxfun(@times,S2t,X2),3) + ...
        sum(bsxfun(@times,S3t(:,:,2:end),X3),3) + sum(bsxfun(@times,S4t,X4),3); %
    
    Srt = sum(bsxfun(@times,S1rt,X1),3) + sum(bsxfun(@times,S2rt,X2),3) + ...
        sum(bsxfun(@times,S3rt(:,:,2:end),X3),3) + sum(bsxfun(@times,S4rt,X4),3); %
    
    clear('X1','X2','X3','X4','X5');
    
    %     Mudança de coordenadas do referêncial(r,t) para o referêncial(x,y),
    %     dos deslocaments e das tensões.
    
    % All pages in T are equal, so the first one is selected to compute the
    % normal cosines
    Ux = cos(T(:,:,1)).*(Ur) - sin(T(:,:,1)).*(Ut);
    Uy = sin(T(:,:,1)).*(Ur) + cos(T(:,:,1)).*(Ut);
    
    % Updating Umax
    Umax = max([Umax,max(max(abs(Ux))),max(max(abs(Uy)))]);
    
    Sx = cos(T(:,:,1)).^2.*Sr + sin(T(:,:,1)).^2.*St - 2.*sin(T(:,:,1)).*cos(T(:,:,1)).*Srt;
    Sy = cos(T(:,:,1)).^2.*St + sin(T(:,:,1)).^2.*Sr + 2.*sin(T(:,:,1)).*cos(T(:,:,1)).*Srt;
    Sxy = sin(T(:,:,1)).*cos(T(:,:,1)).*(Sr-St) + (cos(T(:,:,1)).^2 - sin(T(:,:,1)).^2).*Srt;
    
    % PLOTTING AREA
    
    GlobalX = x(:,:,1) + LocLoop.center(1);
    GlobalY = y(:,:,1) + LocLoop.center(2);
    
    % Writting the fields in a TecPlot compatible file.
    
    % 2 for cycles, one in the x direction, the other in y.
    
    if ~isempty(DirName)
        for jj = 1:NoDiv+1
            for kk = 1:NoDiv+1
                fprintf(FileU,'%0.6e %0.6e %0.6e %0.6e %0.6e %0.6e %0.6e \n',...
                    GlobalX(jj,kk), GlobalY(jj,kk), Ux(jj,kk), Uy(jj,kk),...
                    Sx(jj,kk), Sy(jj,kk), Sxy(jj,kk));
            end
        end
    end
    
    %Ux_Plt;
    subplot(3,2,1);
    contourf(GlobalX, GlobalY, Ux,20,'edgecolor','none'); colormap jet;
    
    %Uy_Plt;
    subplot(3,2,2);
    contourf(GlobalX, GlobalY, Uy,20,'edgecolor','none'); colormap jet;
    
    %Sx_Plt;
    subplot(3,2,3);
    contourf(GlobalX, GlobalY, Sx,20,'edgecolor','none'); colormap jet;
    
    %Sy_Plt;
    subplot(3,2,4);
    contourf(GlobalX, GlobalY, Sy,20,'edgecolor','none'); colormap jet;
    
    %Sxy_Plt;
    subplot(3,2,5);
    contourf(GlobalX, GlobalY, Sxy,20,'edgecolor','none'); colormap jet;
    
end

% While admittedly silly, caxis(caxis) is the workaround that I found for
% it to use a colorbar range taken over all elements rather than over the
% last one alone
subplot(3,2,1); caxis(caxis); colorbar;
subplot(3,2,2); caxis(caxis); colorbar;
subplot(3,2,3); caxis(caxis); colorbar;
subplot(3,2,4); caxis(caxis); colorbar;
subplot(3,2,5); caxis(caxis); colorbar;



% ********************************************************************

% Generating the deformed shape of the structure

% The scaling of the displacements is such that the maximum
% displacement in the structure is represented as 10% of the diagonal
% of the circumscribed rectangle

xmax = (max(Nodes(Loops.nodes(:),1)));
xmin = (min(Nodes(Loops.nodes(:),1)));
ymax = (max(Nodes(Loops.nodes(:),2)));
ymin = (min(Nodes(Loops.nodes(:),2)));

% diagonal of the rectangle circumscribing the structure
D = sqrt((xmax-xmin)^2 + (ymax-ymin)^2);

% the scale factor
S_u = D/(10*Umax);

for ii=1:length(Loops.area)
    
    LocLoop = struct('id',ii,'nodes',Loops.nodes(ii,:),'center',...
        Loops.center(ii,:),'order',Loops.order(ii),'insert',...
        Loops.insert(ii,:),'dim',Loops.dim(ii,:),...
        'materials',Loops.materials(ii,:));
    
    % Static modes x and y
    n = 1:LocLoop.order;
    
    GlobalXn = Nodes(LocLoop.nodes(:),1) ;
    GlobalYn = Nodes(LocLoop.nodes(:),2) ;
    
    
    xn = GlobalXn - LocLoop.center(1);
    yn = GlobalYn - LocLoop.center(2);
    
    % Generating the output grid (local coordinates)
    [xn,yn,N] = ndgrid(xn,yn,n);
    
    % Calculating the R, T coordinates
    Rn = sqrt(xn.^2 + yn.^2);
    Tn = atan2(yn,xn);
    
    % Displacement modes
    U1rn = Rn.^N.*sin((N+1).*Tn);
    U1tn = Rn.^N.*cos((N+1).*Tn);
    U2rn = Rn.^N.*cos((N+1).*Tn);
    U2tn = Rn.^N.*-sin((N+1).*Tn);
    U3rn = Rn.^N.*sin((N-1).*Tn).*((LocLoop.materials(5)*(N-3))+(LocLoop.materials(4)*(N-1)));
    U3tn = Rn.^N.*cos((N-1).*Tn).*((LocLoop.materials(5)*(N+3))+(LocLoop.materials(4)*(N+1)));
    U4rn = Rn.^N.*cos((N-1).*Tn).*((LocLoop.materials(5)*(N-3))+(LocLoop.materials(4)*(N-1)));
    U4tn = -Rn.^N.*sin((N-1).*Tn).*((LocLoop.materials(5)*(N+3))+(LocLoop.materials(4)*(N+1)));
    U1rrbn = 0;
    U1trbn = Rn;
    U2rrbn = cos(Tn);
    U2trbn = -sin(Tn);
    U3rrbn = sin(Tn);
    U3trbn = cos(Tn);
    
    % Vector X que multiplica com os modos (U e S)
    X1(1,1,:) = X(LocLoop.insert(1):LocLoop.insert(1)+LocLoop.dim(1)-1);
    X2(1,1,:) = X(LocLoop.insert(2):LocLoop.insert(2)+LocLoop.dim(2)-1);
    X3(1,1,:) = X(LocLoop.insert(3):LocLoop.insert(3)+LocLoop.dim(3)-1);
    X4(1,1,:) = X(LocLoop.insert(4):LocLoop.insert(4)+LocLoop.dim(4)-1);
    X5 = X(LocLoop.insert(5):LocLoop.insert(5)+LocLoop.dim(5)-1);
    
    % The first page in U3 is cut off from the basis (corresponds to a zero mode)
    Urn = sum(bsxfun(@times,U1rn,X1),3) + sum(bsxfun(@times,U2rn,X2),3) + ...
        sum(bsxfun(@times,U3rn(:,:,2:end),X3),3) + sum(bsxfun(@times,U4rn,X4),3)...
        + U1rrbn(:,:,1)*X5(1) + U2rrbn(:,:,1)*X5(2) + U3rrbn(:,:,1)*X5(3); %
    
    Utn = sum(bsxfun(@times,U1tn,X1),3) + sum(bsxfun(@times,U2tn,X2),3) + ...
        sum(bsxfun(@times,U3tn(:,:,2:end),X3),3) + sum(bsxfun(@times,U4tn,X4),3)...
        + U1trbn(:,:,1)*X5(1) + U2trbn(:,:,1)*X5(2) + U3trbn(:,:,1)*X5(3); %
    
    clear('X1','X2','X3','X4','X5');
    
    % All pages in T are equal, so the first one is selected to compute the
    % normal cosines
    Uxn = cos(Tn(:,:,1)).*(Urn) - sin(Tn(:,:,1)).*(Utn);
    Uyn = sin(Tn(:,:,1)).*(Urn) + cos(Tn(:,:,1)).*(Utn);
    
    
    UxnSc = S_u*diag(Uxn)+GlobalXn;
    UynSc = S_u*diag(Uyn)+GlobalYn;
    
    %Def_Plt;
    subplot(3,2,6);
    def = patch(UxnSc,UynSc,'k'); set(def,'Facecolor','none','Edgecolor','r','LineWidth',2);
    
end

axis tight;
axis off;

if ~isempty(DirName)
    fclose('all');
end

end