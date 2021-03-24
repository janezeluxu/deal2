function [u2]=mapSolution(IBC2,BCval2,p2,IEN_all2,p_All2,vertexData2,p1,IEN_all1,p_All1,vertexData1,variable)
%% Assemble of element level matrix   
Grid_size2 = length(p_All2);
TotalDOF2 = max(max(IEN_all2));
[row2,column2] = K_index(IEN_all2,Grid_size2);

value = zeros(1,length(row2));
rhs = zeros(TotalDOF2,1);
keyCount = 0;

for ele = 1:Grid_size2
%for ele = 1:1
    %% get p2 grid information
    [IENall2,pAll2,xy2] = elementIEN(ele,IEN_all2,p_All2,vertexData2);
    nssl2 = length(IENall2);
    
    n = nIntergerPoints(p2,0);
    qPoints = getintergrate(n);
    
    quadsf = QuadShapeFunc(qPoints,pAll2);
    [ShapeFunc2, ~] = quadsf.QuadShapeFuncTable();
    %% get p1 grid information, shape function evaluate at p2 grid quadrature points
    [IENall1,pAll1,xy1] = elementIEN(ele,IEN_all1,p_All1,vertexData1);
    
    nssl1 = length(IENall1);    
    variable_ele = zeros(nssl1,1);
    for i = 1:nssl1
        variable_ele(i) = variable(IENall1(i));
    end
    %variable_ele

    quadsf = QuadShapeFunc(qPoints,pAll1);
    [ShapeFunc1, ~] = quadsf.QuadShapeFuncTable();
    %% find element Shape Functions, insert stiffness matrix and force vectors
    %  into global matrixs(spars), intergrated at p2 grid quadrature points
    x = xy2(1,:);
    y = xy2(2,:);
    %qPoints
    %ShapeFunc1
    %ShapeFunc2
    %variable_ele
    [elementK,elementF] = eleMapping(nssl2,x,y,qPoints,ShapeFunc1,ShapeFunc2,variable_ele);
    
    %sparse matrix input row, column, and value list, using row and column
    %combination as key
    for i = 1:nssl2
        rowTemp = IENall2(i);
        rhs(rowTemp) = elementF(i) + rhs(rowTemp);
        for j = 1:nssl2           
            %columnTemp = IENall(j);
            %from a given key, find the value index, add to value 
            keyCount = keyCount+1;
            indexKglobal = keyCount;
            %indexKglobal = Kindex(rowTemp,columnTemp);
            value(indexKglobal) = elementK(i,j)+value(indexKglobal);
        end
    end
      
end
%rhs
%value
%[row2,column2,value,rhs] = StrongBC(IBC2,BCval2,row2,column2,value,rhs);
lhs = sparse(row2, column2, value);
%K = full(lhs)
u2 = lhs\rhs;
end

function [elementK,elementF] = eleMapping(sizeN,xCord,yCord,qPoints,ShapeFunc1,ShapeFunc2,u1)
%sizeN = length(IENall2);
elementK = zeros(sizeN, sizeN);
elementF = zeros(sizeN,1);
nP = size(qPoints,1);
for k = 1:nP
    xi = qPoints(k,1);
    eta = qPoints(k,2);
    [JInverse, detJ,gij] = getjacobian(xCord,yCord,xi,eta);
    %detJ
    Jw = detJ*qPoints(k,3);
    uquad = ShapeFunc1(:,k)'*u1;
    elementF = elementF+ShapeFunc2(:,k)*Jw*uquad;
    NbGlobal = ShapeFunc2(:,k);
    
    %diffusion
    elementK = elementK + NbGlobal*NbGlobal'*Jw;
end
%elementK\elementF
end

function [row,column,value,fGlobal] = StrongBC(IBC,BCval,...
    row,column,value,fGlobal)

% Account for BC in forcing vector.

for i = 1:length(column)
    [BCvalIndex,a,BCvalIndex]=intersect(column(i),IBC);
    if BCvalIndex~=0
        fGlobal(row(i))= fGlobal(row(i))-value(i)*BCval(BCvalIndex);
    end
end

%add Boundary conditions from IBC
for i = 1:length(row)
    if (any(row(i) == IBC))
        row(i) = -1;
        column(i) = -1;
        value(i) = NaN;
    end
    if (any(column(i) == IBC))
        row(i) = -1;
        column(i) = -1;
        value(i) = NaN;
    end
end

%add back value as 1 at [IBC(i),IBC(i)]
for i = 1: length (IBC)
    rowTemp = IBC(i);
    fGlobal(rowTemp) = BCval(i);
    %fGlobal(rowTemp) = 0;
    row = [row,rowTemp];
    column = [column,rowTemp];
    value = [value,1];
end


row = row(row~=-1);
column = column(column~=-1);
value = value(~isnan(value));

end