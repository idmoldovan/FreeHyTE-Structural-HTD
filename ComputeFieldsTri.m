function [Ux,Uy,Sx,Sy,Sxy] = ComputeFieldsTri(NoDiv, Nodes, Edges, Loops, X)
% COMPUTEFIELDSTRI computes, stores and plots the final displacement and
% stress fields in the structure.
%
%
% COMPUTEFIELDSTRI is called by MAINTRI. It receives structures Edges and 
% Loops, the solution vector X (calculated in MAINTRI), the node position 
% list Node and the number of divisions for plotting, NoDiv. The fields are
% plotted in a total of (NoDiv+1)^2 points per element. The plotting points 
% belong to the Gauss-Legendre quadrature in every finite element.
% Besides the input arguments, COMPUTEFIELDSTRI loads from STRUCTDEF.mat
% information regarding the file where the results should be stored.
% It returns the arrays Ux and Uy, where the x and y displacement values in 
% the (NoDiv+1)^2 plotting points are stored, and the arrays Sx, Sy and Sxy, 
% that contain the stress values in the same points.
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
% The estimates of the displacement and stress fields are obtained by
% substituting the solution X of the solving system in the domain
% approximations of the respective fields. This is done element by element.
% The estimates are used for:
% * plotting colour maps of all fields and sketching the deformed shape of
% the structure;
% * storing the displacement and stress estimates, calculated in the 
% (NoDiv+1)^2 plotting points, in the Ux, Uy, Sx, Sy and Sxy arrays, which 
% are returned to MAINTRI, for some eventual posterior use;
% * storing the displacement and stress estimates, calculated in the 
% (NoDiv+1)^2 plotting points, in a result file which can be accessed by 
% third party post-processing software (it is pre-formatted for TecPlot). 
% To avoid overloading the disk with (fairly large) undesired result files,
% the result file is only created if the input was saved in the first GUI
% (see Section 5.2.2 of reference [3]).
%
%
% The computation of the field estimates is further covered in Section 7.1
% of reference [2] and Section 4.5 of reference [4]. The output options are
% explained in detail in Section 7.2 of reference [2].

%% Managing the result file
% Loading the save file data from STRUCTDEF.
% * DirName is the folder to save the result file;
% * FileName is the name of the file where the analysis data was saved;
% * EdgesOrder and LoopsOrder are the orders of the approximation bases of
% the edges and elements, as defined in STRUCTDEF GUI.
load('StructDef','DirName','FileName','EdgesOrder','LoopsOrder'); 

% if the user left DirName blank, does not generate written output
if ~isempty(DirName)  
    
    % Generating the file name
    % removing .mat extension
    FileName = FileName(1:end-4);
    % adding information on the refinement orders to the FileName
    FileName = sprintf('%s_ND%dNB%d.dat',FileName,LoopsOrder,EdgesOrder);
    % link the path (DirName) to the file name (FileName). This generates
    % the UFilename string, which allows the creation of the result file
    UFilename = fullfile(DirName,FileName);
    % open result file for writing
    FileU = fopen(UFilename,'w');    
    
    % The following lines are the header needed to prepare the result file
    % to be readable with TecPlot. If you wish to use another software for
    % post-processing, you probably need to change the header.
    % However, the format of the data should be fairly general. The results
    % are written as a matrix of NEL*(NoDiv+1)^2 lines (where NEL is the 
    % total number of finite elements and (NoDiv+1)^2 the total number of 
    % result points per element), and 7 columns. For each result point, the
    % columns list the x and y coordinates, in the global referential, 
    % followed by the values of the displacement and stress fields.
    fprintf(FileU,'TITLE="%s"\n',FileName);
    fprintf(FileU,'VARIABLES="X", "Y", "Ux", "Uy", "Sx", "Sy", "Sxy" \n'); 
    fprintf(FileU,'ZONE T="ND = %d; NB = %d"\n',LoopsOrder,EdgesOrder);

end

% Getting the (global, Cartesian) coordinates of the mesh points
xmesh = [Nodes(Edges.nini(:),1) Nodes(Edges.nfin(:),1)];
ymesh = [Nodes(Edges.nini(:),2) Nodes(Edges.nfin(:),2)];

%% Drawing the mesh, to prepare the plotting of the fields
% The figure Fig contains 6 plots, namely two colour maps of the
% displacement fields in x and y, three colour maps stress fields, and a
% sketch of the deformed shape.
Fig=figure;
set(Fig,'name','Displacement and Stress field','numbertitle','off',...
    'color','w') ;

% Preparing the Ux plot
Ux_Plt = subplot(3,2,1);
hold on; title('X- Displacement');
axis([min(min(xmesh)) max(max(xmesh)) min(min(ymesh)) max(max(ymesh))]); 
% Drawing the mesh, one edge at a time
for ii=1:length(Edges.type)          
    line(xmesh(ii,:),ymesh(ii,:),'color','b');
end
daspect([1 1 1]);

% Preparing the Uy plot
Uy_Plt = subplot(3,2,2);
hold on; title('Y- Displacement');
axis([min(min(xmesh)) max(max(xmesh)) min(min(ymesh)) max(max(ymesh))]);
% Drawing the mesh, one edge at a time
for ii=1:length(Edges.type)
    line(xmesh(ii,:),ymesh(ii,:),'color','b');
end
daspect([1 1 1]);

% Preparing the Sx plot
Sx_Plt = subplot(3,2,3);
hold on; title('X- Stress');
axis([min(min(xmesh)) max(max(xmesh)) min(min(ymesh)) max(max(ymesh))]);
% Drawing the mesh, one edge at a time
for ii=1:length(Edges.type)
    line(xmesh(ii,:),ymesh(ii,:),'color','b');
end
daspect([1 1 1]);

% Preparing the Sy plot
Sy_Plt = subplot(3,2,4);
hold on; title('Y- Stress');
axis([min(min(xmesh)) max(max(xmesh)) min(min(ymesh)) max(max(ymesh))]); 
% Drawing the mesh, one edge at a time
for ii=1:length(Edges.type)
    line(xmesh(ii,:),ymesh(ii,:),'color','b');
end
daspect([1 1 1]);

% Preparing the Sxy plot
Sxy_Plt = subplot(3,2,5);
hold on; title('XY- Stress');
axis([min(min(xmesh)) max(max(xmesh)) min(min(ymesh)) max(max(ymesh))]);
% Drawing the mesh, one edge at a time
for ii=1:length(Edges.type)
    line(xmesh(ii,:),ymesh(ii,:),'color','b');
end
daspect([1 1 1]);

% Preparing the deformed shape sketch
Def_Plt = subplot(3,2,6);
hold on; title('Deformed shape');
% Drawing the mesh, one edge at a time
for ii=1:length(Edges.type)
    line(xmesh(ii,:),ymesh(ii,:),'color',[0.5,0.5,0.5]);
end
daspect([1 1 1]);

% Initializing the maximum displacement (in either direction). It is needed 
% to scale the sketch of the deformed shape of the structure.
Umax = 0;

%% Sweeping the elements to compute the solution fields and draw the colour maps
for ii=1:length(Loops.area) 

    % LocLoop is a local structure where the features of the current
    % element which are directly useful for the calculation of the output
    % fields are stored.
    LocLoop = struct('id',ii,'nodes',Loops.nodes(ii,:),'center',...
        Loops.center(ii,:),'order',Loops.order(ii),'insert',...
        Loops.insert(ii,:),'dim',Loops.dim(ii,:),...
        'materials',Loops.materials(ii,:));

    % Vector containing the orders of the basis
    n(1,1,:) = 1:LocLoop.order;
    
    %% Generating the geometric data
    % Getting coordinates of the nodes of the element (global). 
    LocNodes = Nodes(LocLoop.nodes(:),:);
    
    % Generating the plotting points in the global Cartesian referential.
    % The plotting points belong to the Gauss-Legendre quadrature of the
    % triangular element. 
    [GlobalX,GlobalY,~,~]=triquad(NoDiv+1,LocNodes);
    % Generating the output grid in local coordinates.
    x = GlobalX - LocLoop.center(1);
    y = GlobalY - LocLoop.center(2);
    % Transforming the local Cartesian coordinates into polar.
    r = sqrt(x.^2 + y.^2);  
    t = atan2(y,x);
    % Generating the 3D matrices, for consistency with the programming
    % strategy used in the regular meshes.
    R = repmat(r,[1 1 LocLoop.order]);
    T = repmat(t,[1 1 LocLoop.order]);
    N = repmat(n,[NoDiv+1 NoDiv+1 1]);
    
    %% Computing the basis functions
    % Computing the displacement shape functions, in the polar referential.
    % For a full description of the basis, please refer to Section 3.3.6 of
    % reference [4].
    % U1 basis, generated by odd harmonic potentials
    U1r = R.^N.*sin((N+1).*T);
    U1t = R.^N.*cos((N+1).*T);
    % U2 basis, generated by even harmonic potentials
    U2r = R.^N.*cos((N+1).*T);
    U2t = R.^N.*-sin((N+1).*T);
    % U3 basis, generated by odd biharmonic potentials
    U3r = R.^N.*sin((N-1).*T).*((LocLoop.materials(5)*(N-3))+...
        (LocLoop.materials(4)*(N-1)));
    U3t = R.^N.*cos((N-1).*T).*((LocLoop.materials(5)*(N+3))+...
        (LocLoop.materials(4)*(N+1)));
    % U4 basis, generated by even biharmonic potentials
    U4r = R.^N.*cos((N-1).*T).*((LocLoop.materials(5)*(N-3))+...
        (LocLoop.materials(4)*(N-1)));
    U4t = -R.^N.*sin((N-1).*T).*((LocLoop.materials(5)*(N+3))+...
        (LocLoop.materials(4)*(N+1)));
    % The rigid body modes
    U1rrb = 0;
    U1trb = R;
    U2rrb = cos(T);
    U2trb = -sin(T);
    U3rrb = sin(T);
    U3trb = cos(T);
    
    % Computing the stress shape functions, in the polar referential.
    % For a full description of the basis, please refer to Section 3.3.6 of
    % reference [4].
    % S1 basis, generated by odd harmonic potentials
    S1r =2*LocLoop.materials(5)*N.*R.^(N-1).*sin((N+1).*T);
    S1t = 2*LocLoop.materials(5)*N.*R.^(N-1).*-sin((N+1).*T);
    S1rt = 2*LocLoop.materials(5)*N.*R.^(N-1).*cos((N+1).*T);
    % S2 basis, generated by even harmonic potentials
    S2r = 2*LocLoop.materials(5)*N.*R.^(N-1).*cos((N+1).*T);
    S2t = 2*LocLoop.materials(5)*N.*R.^(N-1).*-cos((N+1).*T);
    S2rt = 2*LocLoop.materials(5)*N.*R.^(N-1).*-sin((N+1).*T);
    % S3 basis, generated by odd biharmonic potentials
    S3r = 2*LocLoop.materials(5)*(LocLoop.materials(4)+...
        LocLoop.materials(5))*N.*R.^(N-1).*((N-3).*sin((N-1).*T));
    S3t = 2*LocLoop.materials(5)*(LocLoop.materials(4)+...
        LocLoop.materials(5))*N.*R.^(N-1).*(-(N+1).*sin((N-1).*T));
    S3rt = 2*LocLoop.materials(5)*(LocLoop.materials(4)+...
        LocLoop.materials(5))*N.*R.^(N-1).*((N-1).*cos((N-1).*T));
    % S4 basis, generated by even biharmonic potentials
    S4r = 2*LocLoop.materials(5)*(LocLoop.materials(4)+...
        LocLoop.materials(5))*N.*R.^(N-1).*((N-3).*cos((N-1).*T));
    S4t = 2*LocLoop.materials(5)*(LocLoop.materials(4)+...
        LocLoop.materials(5))*N.*R.^(N-1).*(-(N+1).*cos((N-1).*T));
    S4rt= 2*LocLoop.materials(5)*(LocLoop.materials(4)+...
        LocLoop.materials(5))*N.*R.^(N-1).*(-(N-1).*sin((N-1).*T));
    
    % Extracting the multipliers of each part of the basis from the X
    % solution vector
    X1(1,1,:) = X(LocLoop.insert(1):LocLoop.insert(1)+LocLoop.dim(1)-1);
    X2(1,1,:) = X(LocLoop.insert(2):LocLoop.insert(2)+LocLoop.dim(2)-1);
    X3(1,1,:) = X(LocLoop.insert(3):LocLoop.insert(3)+LocLoop.dim(3)-1);
    X4(1,1,:) = X(LocLoop.insert(4):LocLoop.insert(4)+LocLoop.dim(4)-1);
    X5 = X(LocLoop.insert(5):LocLoop.insert(5)+LocLoop.dim(5)-1);
    
    %% Computing the displacement and stress fields
    % Computing the displacement and stres fields in the polar referential.
    % They are the product of the basis functions with the corresponding
    % solution vectors. The first function in U3&S3 is cut off from the 
    % basis (corresponds to a rigid mode).
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
    
    % Clearing the Xi variables for reuse in the next element
    clear('X1','X2','X3','X4','X5','n');
    
    % Transforming the displacement field from the polar to the Cartesian
    % referential.
    % All pages in T are equal, so the first one is selected to compute the
    % normal cosines.
    Ux = cos(T(:,:,1)).*(Ur) - sin(T(:,:,1)).*(Ut);
    Uy = sin(T(:,:,1)).*(Ur) + cos(T(:,:,1)).*(Ut);
    
    % Updating the maximum structural displacement, Umax
    Umax = max([Umax,max(max(abs(Ux))),max(max(abs(Uy)))]);
    
    % Transforming the stress field from the polar to the Cartesian
    % referential.
    Sx = cos(T(:,:,1)).^2.*Sr + sin(T(:,:,1)).^2.*St - ...
        2.*sin(T(:,:,1)).*cos(T(:,:,1)).*Srt;
    Sy = cos(T(:,:,1)).^2.*St + sin(T(:,:,1)).^2.*Sr + ...
        2.*sin(T(:,:,1)).*cos(T(:,:,1)).*Srt;
    Sxy = sin(T(:,:,1)).*cos(T(:,:,1)).*(Sr-St) + ...
        (cos(T(:,:,1)).^2 - sin(T(:,:,1)).^2).*Srt;
    
    %% Storing the displacement and stress fields in the result file 
    % Writing the fields in the TecPlot compatible file, if requested by
    % the user. The results are written as a matrix of NEL*(NoDiv+1)^2 lines 
    % (where NEL is the total number of finite elements and (NoDiv+1)^2 the 
    % total number of result points per element), and 7 columns. For each 
    % result point, the columns list the x and y coordinates, in the global 
    % referential, followed by the values of the displacement and stress 
    % fields.
    if ~isempty(DirName)
        for jj = 1:NoDiv+1
            for kk = 1:NoDiv+1
                fprintf(FileU,'%0.6e %0.6e %0.6e %0.6e %0.6e %0.6e %0.6e \n',...
                    GlobalX(jj,kk), GlobalY(jj,kk), Ux(jj,kk), Uy(jj,kk),...
                    Sx(jj,kk), Sy(jj,kk), Sxy(jj,kk));
            end
        end
    end
    
    %% Plotting the displacement and stress fields as colour maps
    %Ux_Plt;
    subplot(3,2,1);
    contourf(GlobalX, GlobalY, Ux,20,'edgecolor','none'); colormap jet;
    
    %Uy_Plt;
    subplot(3,2,2);
    contourf(GlobalX, GlobalY, Uy,20,'edgecolor','none'); colormap jet;
    
    %Sx_Plt;
    subplot(3,2,3);
    contourf(GlobalX, GlobalY, Sx,20,'edgecolor','none');  colormap jet;
    
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
subplot(3,2,3); caxis(caxis); colorbar; % caxis([-300,300])
subplot(3,2,4); caxis(caxis); colorbar;
subplot(3,2,5); caxis(caxis); colorbar;


%% ********************************************************************

%% Generating the deformed shape of the structure
% The deformed shape of the structure is sketched by connecting with
% straight lines the displacements of the vertices of the elements. In
% order to obtain the displacements of the vertices, the same procedure
% used for the drawing of the color maps is applied, this time on the nodes
% of the elements only. It is noted that the colour maps do not include the
% vertices of the elements because they do not belong to the Gauss-Legendre
% quadrature.

% The displacements are scaled such that the maximum displacement in the 
% structure is represented as roughly 10% of the diagonal of the rectangle
% that circumscribes the structure.
% Getting the vertices of the rectangle that circumscribes the structure.
xmax = (max(Nodes(Loops.nodes(:),1)));
xmin = (min(Nodes(Loops.nodes(:),1)));
ymax = (max(Nodes(Loops.nodes(:),2)));
ymin = (min(Nodes(Loops.nodes(:),2)));
% Computing the diagonal of the rectangle circumscribing the structure
D = sqrt((xmax-xmin)^2 + (ymax-ymin)^2);
% Computing the scale factor
S_u = D/(10*Umax);

% Sweeping the elements to compute their nodal displacements
for ii=1:length(Loops.area)
    
    % LocLoop is a local structure where the features of the current
    % element which are directly useful for the calculation of the output
    % fields are stored.
    LocLoop = struct('id',ii,'nodes',Loops.nodes(ii,:),'center',...
        Loops.center(ii,:),'order',Loops.order(ii),'insert',...
        Loops.insert(ii,:),'dim',Loops.dim(ii,:),...
        'materials',Loops.materials(ii,:));
    
    % Vector containing the orders of the basis
    n = 1:LocLoop.order;
    
    % Getting the global coordinates of the nodes...
    GlobalXn = Nodes(LocLoop.nodes(:),1) ;
    GlobalYn = Nodes(LocLoop.nodes(:),2) ;
    % ...and transforming them to the local referential      
    xn = GlobalXn - LocLoop.center(1);
    yn = GlobalYn - LocLoop.center(2);
    % Generating the output grid (local coordinates)
    [xn,yn,N] = ndgrid(xn,yn,n);
    
    % Transforming the local Cartesian coordinates into polar.
    Rn = sqrt(xn.^2 + yn.^2);
    Tn = atan2(yn,xn);
    
    %% Computing the basis functions
    % Computing the displacement shape functions, in the polar referential.
    % For a full description of the basis, please refer to Section 3.3.6 of
    % reference [4].
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
    
    % Extracting the multipliers of each part of the basis from the X
    % solution vector
    X1(1,1,:) = X(LocLoop.insert(1):LocLoop.insert(1)+LocLoop.dim(1)-1);
    X2(1,1,:) = X(LocLoop.insert(2):LocLoop.insert(2)+LocLoop.dim(2)-1);
    X3(1,1,:) = X(LocLoop.insert(3):LocLoop.insert(3)+LocLoop.dim(3)-1);
    X4(1,1,:) = X(LocLoop.insert(4):LocLoop.insert(4)+LocLoop.dim(4)-1);
    X5 = X(LocLoop.insert(5):LocLoop.insert(5)+LocLoop.dim(5)-1);
    
    %% Computing the displacement and stress fields
    % Computing the displacement and stres fields in the polar referential.
    % They are the product of the basis functions with the corresponding
    % solution vectors. The first function in U3 is cut off from the
    % basis (corresponds to a rigid mode).
    Urn = sum(bsxfun(@times,U1rn,X1),3) + sum(bsxfun(@times,U2rn,X2),3) + ...
        sum(bsxfun(@times,U3rn(:,:,2:end),X3),3) + sum(bsxfun(@times,U4rn,X4),3)...
        + U1rrbn(:,:,1)*X5(1) + U2rrbn(:,:,1)*X5(2) + U3rrbn(:,:,1)*X5(3); 
    
    Utn = sum(bsxfun(@times,U1tn,X1),3) + sum(bsxfun(@times,U2tn,X2),3) + ...
        sum(bsxfun(@times,U3tn(:,:,2:end),X3),3) + sum(bsxfun(@times,U4tn,X4),3)...
        + U1trbn(:,:,1)*X5(1) + U2trbn(:,:,1)*X5(2) + U3trbn(:,:,1)*X5(3); 
    
    % Clearing the Xi variables for reuse in the next element
    clear('X1','X2','X3','X4','X5');
    
    % Transforming the displacement field from the polar to the Cartesian
    % referential.
    % All pages in T are equal, so the first one is selected to compute the
    % normal cosines.
    Uxn = cos(Tn(:,:,1)).*(Urn) - sin(Tn(:,:,1)).*(Utn);
    Uyn = sin(Tn(:,:,1)).*(Urn) + cos(Tn(:,:,1)).*(Utn);
    
    % Scaling the nodal displacements using the scaling factor S_u
    UxnSc = S_u*diag(Uxn)+GlobalXn;
    UynSc = S_u*diag(Uyn)+GlobalYn;
    
    %% Sketching the deformed shape of the strucutre
    subplot(3,2,6);
    def = patch(UxnSc,UynSc,'k'); set(def,'Facecolor','none','Edgecolor','r','LineWidth',2);

end
axis tight;
axis off;

%% House cleaning
if ~isempty(DirName)
    fclose('all');
end
end