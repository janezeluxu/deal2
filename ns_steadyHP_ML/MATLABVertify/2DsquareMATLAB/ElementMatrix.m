function [elementL,elementR] = ElementMatrix(xcord,ycord,va,miu,nInt)

%miu = 1;
%nInt = 3;

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
            
Ruconv = zeros(4,1);
Rudiff = zeros(4,1);
Rumtauc = zeros(4,1);
Rvconv = zeros(4,1);
Rvdiff = zeros(4,1);
Rctaum = zeros(4,1);
Rcc = zeros(4,1);
Rutaum = zeros(4,1);
Rvtaum = zeros(4,1);
Rvmtauc = zeros(4,1);

nsf = 4;
K11 = zeros(nsf, nsf);
K12 = zeros(nsf, nsf);
K21 = zeros(nsf, nsf);
K22 = zeros(nsf, nsf);
G1 = zeros(nsf, nsf);
G2 = zeros(nsf, nsf);
C = zeros(nsf, nsf);
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
   ubar = uele-tauM*rum;
   vbar = vele-tauM*rvm;
   rc = ugrad(1)+vgrad(2);
   
   gradNx = gradNaGlobal(:,1);
   gradNy = gradNaGlobal(:,2);
     
   Ruconv = Ruconv+ShapeFunc(:,k)*(ubar*ugrad(1)+vbar*ugrad(2))*Jw(k)-gradNx*pele*Jw(k);
   Rudiff= Rudiff+miu*gradNaGlobal*[ugrad(1);ugrad(2)]*Jw(k);
   Rutaum = Rutaum+(rum*tauM*[uele,vele]*gradNaGlobal')'*Jw(k);
   Rumtauc = Rumtauc+rc*tauC*gradNx*Jw(k);
   
   Rvconv = Rvconv+ShapeFunc(:,k)*(ubar*vgrad(1)+vbar*vgrad(2))*Jw(k)-gradNy*pele*Jw(k);
   Rvdiff= Rvdiff+miu*gradNaGlobal*[vgrad(1);vgrad(2)]*Jw(k);
   Rvtaum = Rvtaum+(rvm*tauM*[uele,vele]*gradNaGlobal')'*Jw(k);
   Rvmtauc = Rvmtauc+rc*tauC*gradNy*Jw(k);
   
   Rcc = Rcc+ShapeFunc(:,k)*(ugrad(1)+vgrad(2))*Jw(k);
   Rctaum = Rctaum+tauM*(gradNx*rum+gradNy*rvm)*Jw(k);
   % element matrix component 1-13
   SF = ShapeFunc(:,k);
   gradNGlobal = gradNaGlobal';
   rou = 1;
   K11 = K11+SF*rou*([ubar,vbar]*gradNGlobal)*Jw(k);
   K11 = K11+gradNx*miu*gradNx'*Jw(k)+gradNy*miu*gradNy'*Jw(k);
   K11 = K11+gradNGlobal'*tauM*[uele;vele]*([uele,vele]*gradNGlobal)*Jw(k);
   K11 = K11+(gradNx*tauC*gradNx')*Jw(k);
   
   K12 = K12+(gradNx*tauC*gradNy')*Jw(k);%miu*gradNy*(gradNx')*Jw(k);
   
   K21 = K21+(gradNy*tauC*gradNx')*Jw(k);%+miu*gradNx*(gradNy')*Jw(k);
   
   K22 = K22+SF*rou*([ubar,vbar]*gradNGlobal)*Jw(k);
   K22 = K22+gradNx*miu*gradNx'*Jw(k)+gradNy*miu*gradNy'*Jw(k);
   K22 = K22+gradNGlobal'*tauM*[uele;vele]*([uele,vele]*gradNGlobal)*Jw(k);
   K22 = K22+(gradNy*tauC*gradNy')*Jw(k);
   
   G1 = G1-(gradNx*SF')*Jw(k);
   G2 = G2-(gradNy*SF')*Jw(k);
   
    
   C = C+(gradNx*tauM*gradNx')*Jw(k)+(gradNy*tauM*gradNy')*Jw(k);
end

Rum = Ruconv+Rudiff+Rutaum+Rumtauc;
Rvm = Rvconv+Rvdiff+Rvtaum+Rvmtauc;
Rc = Rcc+Rctaum;
D1 = -G1';
D2 = -G2';
elementL = [K11,K12,G1;K21,K22,G2;D1,D2,C];
%Rc
elementR = -[Rum;Rvm;Rc];

%Rm = Rum^2+Rvm^2;
eL = [K11(1,1),K12(1,1),G1(1,1),K11(1,2),K12(1,2),G1(1,2),K11(1,4),K12(1,4),G1(1,4),K11(1,3),K12(1,3),G1(1,3);
    K21(1,1),K22(1,1),G2(1,1),K21(1,2),K22(1,2),G2(1,2),K21(1,4),K22(1,4),G2(1,4),K21(1,3),K22(1,3),G2(1,3);
    D1(1,1),D2(1,1),C(1,1),D1(1,2),D2(1,2),C(1,2),D1(1,4),D2(1,4),C(1,4),D1(1,3),D2(1,3),C(1,3);
    K11(2,1),K12(2,1),G1(2,1),K11(2,2),K12(2,2),G1(2,2),K11(2,4),K12(2,4),G1(2,4),K11(2,3),K12(2,3),G1(2,3);
    K21(2,1),K22(2,1),G2(2,1),K21(2,2),K22(2,2),G2(2,2),K21(2,4),K22(2,4),G2(2,4),K21(2,3),K22(2,3),G2(2,3);
    D1(2,1),D2(2,1),C(2,1),D1(2,2),D2(2,2),C(2,2),D1(2,4),D2(2,4),C(2,4),D1(2,3),D2(2,3),C(2,3);
    K11(4,1),K12(4,1),G1(4,1),K11(4,2),K12(4,2),G1(4,2),K11(4,4),K12(4,4),G1(4,4),K11(4,3),K12(4,3),G1(4,3);
    K21(4,1),K22(4,1),G2(4,1),K21(4,2),K22(4,2),G2(4,2),K21(4,4),K22(4,4),G2(4,4),K21(4,3),K22(4,3),G2(4,3);
    D1(4,1),D2(4,1),C(4,1),D1(4,2),D2(4,2),C(4,2),D1(4,4),D2(4,4),C(4,4),D1(4,3),D2(4,3),C(4,3)
    K11(3,1),K12(3,1),G1(3,1),K11(3,2),K12(3,2),G1(3,2),K11(3,4),K12(3,4),G1(3,4),K11(3,3),K12(3,3),G1(3,3);
    K21(3,1),K22(3,1),G2(3,1),K21(3,2),K22(3,2),G2(3,2),K21(3,4),K22(3,4),G2(3,4),K21(3,3),K22(3,3),G2(3,3);
    D1(3,1),D2(3,1),C(3,1),D1(3,2),D2(3,2),C(3,2),D1(3,4),D2(3,4),C(3,4),D1(3,3),D2(3,3),C(3,3);];
end