classdef SquareShapeFunc < handle
%calculate shape functions for Bernstain polynomial in lambda space,
%simplex order can be uniform or variable/transition elements
properties (Access = private)
        %variable p element IEN
        elementIEN;        
        %uniform element Order
        Order;             
        %variable p element order list
        OrderAll;          
        %list of quadrature points
        QuadraturePoints;

    end

methods(Access = public)
    function sqaure = SquareShapeFunc(QPoints,varargin)
        sqaure.QuadraturePoints = QPoints;
        %if the simplex has uniform order, only need p information
        if length(varargin) ==1
            sqaure.Order = varargin{1};
            sqaure.elementIEN = nan;
            sqaure.OrderAll = nan;
        %if element is variable order, will need the element IEN array and order list   
        else
            sqaure.Order = nan;
            sqaure.elementIEN = varargin{1};
            sqaure.OrderAll = varargin{2};
        end
    end
    
    function [ShapeFuncTable, divShapeFuncTable] = uniformPShapeFunctionTable(square)
        %this function will call uniformPShapeFunction and give the shape
        %function table of a quadrature list 
        qPoints = square.QuadraturePoints;
        p = square.Order;
        nShapeFunc = 4;
        nQuadPoints = size(qPoints,1);
        ShapeFuncTable = zeros(nShapeFunc,nQuadPoints);
        divShapeFuncTable = zeros(nShapeFunc,2,nQuadPoints);
        for i = 1: nQuadPoints
            IPoint = qPoints(i,:);
            [ ShapeFunc, divShapeFunc ] = uniformPShapeFunction(square,IPoint);
            ShapeFuncTable(:,i) =  ShapeFunc;
            divShapeFuncTable(:,:,i) = divShapeFunc;
        end
        
    end
    
    function [ ShapeFunc, divShapeFunc ] = uniformPShapeFunction(square,IntPoint)
        %calculate the uniform order shape functions at one given intergral point, this
        %function calls ShapeFunction
        p = square.Order;
        
        [ ShapeFunc, divShapeFunc ] = square.ShapeFunction(IntPoint);
    end
    
end

methods(Static)
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
    
end
end