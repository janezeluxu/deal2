function getvalue()
v1 = [0.9875 0.9875];
 v2 = [1 0.9875];
 v3 = [0.9875 1];
 v4 = [1 1 ];

vPoint =[ 0.994365 0.994365];

IntPoint = mapToInt(vPoint,v1,v2,v3,v4);

[ ShapeFunc, divShapeFunc ] = ShapeFunction(IntPoint);
ShapeFunc
u = [0.25409,0,1,1];
uPoint = u*ShapeFunc
end


function IntPoint = mapToInt(vPoint,v1,v2,v3,v4)
xi = (vPoint(1)-v1(1))/(v2(1)-v1(1))
eta = (vPoint(2)-v1(2))/(v3(2)-v1(2))
IntPoint(1) = xi;
 IntPoint(2) = eta;
end

function [ ShapeFunc, divShapeFunc ] = ShapeFunction(IntPoint)
        %give one shape function and derivative value at one intergral point
        %and one lexico order combination in lambda space
        xi = IntPoint(1);
        eta = IntPoint(2);
        
        ShapeFunc = [(1-xi)*(1-eta);xi*(1-eta);xi*eta;(1-xi)*eta];
        divShapeFunc = [-(1-eta),-(1-xi);1-eta,-xi;eta,xi;-eta,(1-xi)];
        
        %if (a1+a2+a3 == 1)
        %    divShapeFunc = [-1,-1;1,0;0,1];
        %end
        
    end