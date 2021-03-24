function [L2_ele,H1_ele] = estimateError(variable,meshData,IEN,vertexData)
global miu;
global Grid_size;

L2_ele = zeros(Grid_size,1);
H1_ele = zeros(Grid_size,1);
nInt = 3;
    [xi, w] = GaussQuad(nInt,1);
lxi = length(xi);
quadraturePoints = zeros(lxi*lxi,3);
for i = 1:lxi
    for j = 1:lxi
        n = (i-1)*lxi+j;
        quadraturePoints(n,:)=[xi(i),xi(j),w(i)*w(j)];
    end
end
%L2_u = 0;
c = 2/(3^0.5);
for ele = 1:Grid_size 
%for ele = 1:1 
    ele;
    [IENall,xcord,ycord] = elementIEN(ele,meshData,vertexData,IEN);
    nssl = length(IENall);
    va = zeros(nssl,1);
    for i = 1:nssl
        va(i) = variable(IENall(i));
    end
    nP = size(quadraturePoints,1);
    %quadraturePoints
    [ShapeFuncTable,divSFtable] = ShapeTable(nInt,1);
    ShapeFunc = ShapeFuncTable{1};
    divShapeFunc = divSFtable{1};
    J = [(xcord(2)-xcord(1)),0;0,(ycord(4)-ycord(1))];
    JInverse = J^-1;
    gij = [(JInverse(1,1))^2,0;0,(JInverse(2,2))^2];
    %gij = 4*gij;
    detJ = det(J);
    Jw = detJ*quadraturePoints(:,3);
    
    nsf = size(ShapeFunc,1);
    u = va(1:nsf);
    v = va(nsf+1:2*nsf);
    p = va(2*nsf+1:3*nsf);
    L2_u = 0;
    H1_u = 0;
    for k = 1:nP
        uele = u'*ShapeFunc(:,k);
        vele = v'*ShapeFunc(:,k);
        pele = p'*ShapeFunc(:,k);
        
        C2 = 9;
        tau1sqinv = [uele,vele]*gij*[uele;vele];
        A = gij'*gij;
        tau2sqinv = C2*miu^2*(A(1,1)+A(2,2));
        tauM = 1/sqrt(tau1sqinv+tau2sqinv);
        tauC = (1)/(tauM*(gij(1,1)+gij(2,2)));
        
        %ShapeFunc(:,k)
        
        gradNaGlobal = divShapeFunc(:,:,k)*JInverse;
        ugrad = u'*gradNaGlobal;
        vgrad = v'*gradNaGlobal;
        pgrad = p'*gradNaGlobal;
        
        rum = uele*ugrad(1)+vele*ugrad(2)+pgrad(1);
        rvm = uele*vgrad(1)+vele*vgrad(2)+pgrad(2);
        magrm = rum^2+rvm^2;
        %tauM
        %detJ
        H1c = miu^-0.5*tauM^0.5;
        H1_u = H1_u + H1c*H1c*magrm*quadraturePoints(k,3)*detJ;
        L2_u = L2_u + c*tauM*c*tauM*magrm*quadraturePoints(k,3)*detJ;
    end

    L2_ele(ele) = L2_u^0.5;
    H1_ele(ele) = H1_u^0.5;
end
%L2_u = L2_u^0.5
end