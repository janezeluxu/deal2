function [L2_ele,H1_ele]=getErrorwithRef(variable,refSolutionfile,Nsol,Msol)
global miu;
global Grid_size;
nInt = 3;
[ShapeFuncTable,divSFtable] = ShapeTable(nInt,1);
ShapeFunc = ShapeFuncTable{1};
divShapeFunc = divSFtable{1};

N = 81;
M = 81;
[meshDataRef,IENRef,vertexDataRef,IBCRef,BCvalRef,TotalNodeRef,totalDOFRef] = ...
    createMesh(N,M);
totalEleRef = (N-1)*(M-1);

variableRef = load(refSolutionfile);
%N = 80;
%M = 80;
[meshData,IEN,vertexData,IBC,BCval,TotalNode,totalDOF] = createMesh(Nsol,Msol);
Error_u = 0;
Error_v = 0;
Error_p = 0;
Error_ux = 0;
Error_uy = 0;
L2_ele = zeros(Grid_size,1);
H1_ele = zeros(Grid_size,1);

for ele = 1:Grid_size 
%for ele = 1:1 
    %ele
    [IENall,xcord,ycord] = elementIEN(ele,meshData,vertexData,IEN);
    nssl = length(IENall);
    va = zeros(nssl,1);
    for i = 1:nssl
        va(i) = variable(IENall(i));
    end
    
    nsf = size(ShapeFunc,1);
    u = va(1:nsf);
    v = va(nsf+1:2*nsf);
    p = va(2*nsf+1:3*nsf);
    
    J = [(xcord(2)-xcord(1)),0;0,(ycord(4)-ycord(1))];
    JInverse = J^-1;
    
    [xi, w] = GaussQuad(nInt,1);
    lxi = length(xi);
    quadraturePoints = zeros(lxi*lxi,3);
    for i = 1:lxi
        for j = 1:lxi
            n = (i-1)*lxi+j;
            quadraturePoints(n,:)=[xi(i),xi(j),w(i)*w(j)];
        end
    end
    nP = size(quadraturePoints,1);
    
    detJ = det(J);
    Jw = detJ*quadraturePoints(:,3);
    Error_u = 0;
    Error_v = 0;
    for k = 1:nP
        
        uele = u'*ShapeFunc(:,k);
        vele = v'*ShapeFunc(:,k);
        pele = p'*ShapeFunc(:,k);
        
        gradNaGlobal = divShapeFunc(:,:,k)*JInverse;
        ugrad = u'*gradNaGlobal;
        vgrad = v'*gradNaGlobal;
        pgrad = p'*gradNaGlobal;
   
        %% map to reference solution
        %xcord
        %ycord
        [x,y] = getrealcoordinate(xcord,ycord,quadraturePoints(k,:));
        eleRef = findRefelement(x,y,meshDataRef,IENRef,vertexDataRef,totalEleRef);
        [IENallref,xcordRef,ycordRef] = elementIEN(eleRef,meshDataRef,vertexDataRef,IENRef);
        %xcordRef
        %ycordRef
        evalPoint = maptoparametric(x,y,xcordRef,ycordRef);
        
        Jref = [(xcordRef(2)-xcordRef(1)),0;0,(ycordRef(4)-ycordRef(1))];
        JInverseref = Jref^-1;
    
        nssl = length(IENallref);
        vaRef = zeros(nssl,1);
        for i = 1:nssl
            vaRef(i) = variableRef(IENallref(i));
        end
        
        [ ShapeFuncRef, divShapeFuncRef ] = getShapeFunc(evalPoint);
        nsf = size(ShapeFuncRef,1);
        uRef = vaRef(1:nsf);
        vRef = vaRef(nsf+1:2*nsf);
        pRef = vaRef(2*nsf+1:3*nsf);
        
        ueleref = uRef'*ShapeFuncRef;
        veleref = vRef'*ShapeFuncRef;
        peleref = pRef'*ShapeFuncRef;
        
        gradNaGlobal = divShapeFuncRef*JInverseref;
        ugrad_ref = u'*gradNaGlobal;
        vgrad_ref = v'*gradNaGlobal;
        pgrad_ref = p'*gradNaGlobal;
   
        Error_ux = (((ugrad(1) - ugrad_ref(1))).^2*quadraturePoints(k,3)*detJ)+Error_ux;
        Error_uy = (((ugrad(2) - ugrad_ref(2))).^2*quadraturePoints(k,3)*detJ)+Error_uy;
        
        %ueleref
        %veleref
        
        %uele
        %vele
        %pele
        
        %peleref
        %% get true error
        Error_u = (((uele - ueleref)).^2*quadraturePoints(k,3)*detJ)+Error_u;
        Error_v = (((vele - veleref)).^2*quadraturePoints(k,3)*detJ)+Error_v;
        Error_p = (((pele - peleref)).^2*quadraturePoints(k,3)*detJ)+Error_p;
    end
    %Error_u
    %Error_v
    H1_ele(ele) = (Error_ux+Error_uy)^0.5;
    L2_ele(ele) = (Error_u+Error_v)^0.5;
end

% Error_uxy = (Error_ux+Error_uy)^0.5
% Error_u = Error_u^0.5
% Error_v = Error_v^0.5;
% Error_p = Error_p^0.5;
end

function [ ShapeFunc, divShapeFunc ] = getShapeFunc(IntPoint)
%give one shape function and derivative value at one intergral point
%and one lexico order combination in lambda space
xi = IntPoint(1);
eta = IntPoint(2);

ShapeFunc = [(1-xi)*(1-eta);xi*(1-eta);xi*eta;(1-xi)*eta];
divShapeFunc = [-(1-eta),-(1-xi);1-eta,-xi;eta,xi;-eta,(1-xi)];


end

function [xc,yc] = getrealcoordinate(x,y,IntPoint)

[ ShapeFunc ] = getShapeFunc(IntPoint);
N1 = ShapeFunc(1);
N2 = ShapeFunc(2);
N3 = ShapeFunc(3);
N4 = ShapeFunc(4);
xc = N1*x(1)+N2*x(2)+N3*x(3)+N4*x(4);
yc = N1*y(1)+N2*y(2)+N3*y(3)+N4*y(4);
end

function eleRef = findRefelement(x,y,meshDataRef,IENRef,vertexDataRef,totalEleRef)
distanceToCenter = zeros(totalEleRef,1);
for ele = 1:totalEleRef
    [IENall,xcord,ycord] = elementIEN(ele,meshDataRef,vertexDataRef,IENRef);
    x1 = xcord(1);
    x2 = xcord(2);
    x3 = xcord(3);
    x4 = xcord(4);
    
    y1 = ycord(1);
    y2 = ycord(2);
    y3 = ycord(3);
    y4 = ycord(4);
    
    xmid = (x1+x2+x3+x4)/4;
    ymid = (y1+y2+y3+y4)/4;
    
    distanceToCenter(ele) = sqrt((x-xmid)^2+(y-ymid)^2);
end
[M,eleRef]=min(distanceToCenter);
end

function EvalPoints = maptoparametric(x,y,xCordM1,yCordM1)
x1 = xCordM1(1);
x2 = xCordM1(2);
x3 = xCordM1(3);
x4 = xCordM1(4);

y1 = yCordM1(1);
y2 = yCordM1(2);
y3 = yCordM1(3);
y4 = yCordM1(4);

%% get evalPoints in xi,eta coordinate

a = (y4-y1)/(x4-x1);
b = y1-a*x1;
l1 = abs(a*x-y+b)/sqrt(a^2+1);
l = sqrt((x1-x2)^2+(y2-y1)^2);
xi = l1/l;

a = (y2-y1)/(x2-x1);
b = y1-a*x1;
m1 = abs(a*x-y+b)/sqrt(a^2+1);
m = sqrt((x1-x4)^2+(y4-y1)^2);
eta = m1/m;

xi = (x-x1)/(x2-x1);
eta = (y-y1)/(y4-y1);
EvalPoints(1) = abs(xi);
EvalPoints(2) = abs(eta);
end