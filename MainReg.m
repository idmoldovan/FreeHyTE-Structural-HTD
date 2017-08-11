function MainReg
close all;

% Main processor program

% Launch pre-processing
[NGP, Nodes, Edges, Loops, BConds, NoDiv] = InputProcReg;

% Assign parts to elements and sides in the global system
[Edges,Loops,Dim] = AssignParts(Edges, Loops, BConds);
LHS = zeros(Dim);
RHS = zeros(Dim,1);

% Initialization of Gauss weights & abscissas (on a -1:1 interval)
[abscissa,weight] = gauleg(NGP, -1, 1);

% Generating & allocating the LHS matrices
LHS = K11(Edges, Loops, LHS, abscissa, weight);
LHS = K12(Edges, Loops, LHS, abscissa, weight);
LHS = K13(Edges, Loops, LHS, abscissa, weight);
LHS = K14(Edges, Loops, LHS, abscissa, weight);
LHS = K22(Edges, Loops, LHS, abscissa, weight);
LHS = K23(Edges, Loops, LHS, abscissa, weight);
LHS = K24(Edges, Loops, LHS, abscissa, weight);
LHS = K33(Edges, Loops, LHS, abscissa, weight);
LHS = K34(Edges, Loops, LHS, abscissa, weight);
LHS = K44(Edges, Loops, LHS, abscissa, weight);
LHS = B1(Edges, Loops, LHS, abscissa, weight);
LHS = B2(Edges, Loops, LHS, abscissa, weight);
LHS = B3(Edges, Loops, LHS, abscissa, weight);
LHS = B4(Edges, Loops, LHS, abscissa, weight);
LHS = B5(Edges, Loops, LHS, abscissa, weight);

% Generating & allocating the RHS vectors
RHS = T1(Edges, Loops, BConds, RHS, abscissa, weight);
RHS = T2(Edges, Loops, BConds, RHS, abscissa, weight);
RHS = T3(Edges, Loops, BConds, RHS, abscissa, weight);
RHS = T4(Edges, Loops, BConds, RHS, abscissa, weight);
RHS = T5(Edges, Loops, BConds, RHS, abscissa, weight);
RHS = U(Edges, BConds, RHS, abscissa, weight);

% Scaling the system
Sc = sqrt(diag(LHS)).^-1;
Sc(isinf(Sc)) = 1;
Sc = diag(Sc);

SparseQ = 0;
try % checking if out of memory...
    ScLHS = Sc' * LHS * Sc;
catch % ... and going sparse if it is
    Sc = sparse(Sc);
    LHS = sparse(LHS);
    ScLHS = Sc' * LHS * Sc;
    SparseQ = 1;
end

ScRHS = Sc' * RHS;

% Solving the FE system
if ~SparseQ % if not sparse
    CndNo=rcond(ScLHS);
    if (CndNo<eps)
        warning('local:NumericalChk',...
            'Ststem condition number is %0.5g.\n',...
            CndNo);
        ScX = pinv(ScLHS)*ScRHS;
    else
        ScX = ScLHS\ScRHS;
    end
else % if sparse
    ScX = ScLHS\ScRHS;
end

% Scaling back
X = Sc * ScX;

% Constructing the displacement field 
[Ux, Uy, Sx, Sy, Sxy] = ComputeFieldsReg(NoDiv, Nodes, Edges, Loops,X);

% restoring warnings
warning('on','MATLAB:DELETE:FileNotFound');
warning('on','MATLAB:load:variableNotFound');

end