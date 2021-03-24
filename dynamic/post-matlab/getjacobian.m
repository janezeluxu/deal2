function [JInverse, detJ,gij] = getjacobian(xCord,yCord,xi,eta)

J = [0,0;0,0];
dNl1 = [-(1-eta), 1-eta, -eta,eta]; % derivative of N1, N2 and N4 with xi
dNl2 = [-(1-xi), -xi, (1-xi),xi]; % derivative of N1, N2 and N4 with eta

%divShapeFunc = [-(1-eta),-(1-xi);1-eta,-xi;eta,xi;-eta,(1-xi)];
for i = 1:4
    J(1,1) = J(1,1)+xCord(i)*dNl1(i); % derivative of x with xi
    J(1,2) = J(1,2)+xCord(i)*dNl2(i); % derivative of x with eta
    J(2,1) = J(2,1)+yCord(i)*dNl1(i); % derivative of y with xi
    J(2,2) = J(2,2)+ yCord(i)*dNl2(i); % derivative of y with eta
end
JInverse = J^-1;
gij = JInverse'*JInverse;
detJ = det(J);
end