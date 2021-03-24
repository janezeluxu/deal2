function [concentration]=LinePlot(mesh_file,u,x,y,whichvariable)
 global Grid_size;

[solution_ien,order_list,vertexData] = readin_mesh(mesh_file);
Grid_size = size(solution_ien,1);
 %u = load(solution);
 
 concentration = zeros(1,length(y));
 x1 = x(1)
 y1 = y(1)
 %% start the search from boundary elements
 [startingEle]=findNextEle(x(1),y(1),1:Grid_size,...
        solution_ien,order_list,vertexData);
 %startingEle = 12461;
 lastEle = startingEle;
 
 eleList = 1:Grid_size;%[35,36,37,38,39,451,452,453,454,455,541,542,543,...
    % 548,549,550,558,559,561,570,572,574]+1;
 for i = 1:length(x)
     %lastEle
     xx = x(i);
     yy = y(i);
    %eleSearch = getSurrondEle(lastEle,meshData,...
     %   solution_ien,order_list,vertexData);
    [newEle] = findNextEle(x(i),y(i),eleList,...
        solution_ien,order_list,vertexData);
    %x(i),y(i)
    concentration(i) = getSolution(u,x(i),y(i),newEle,...
        solution_ien,order_list,vertexData,whichvariable);
    %cc = concentration(i)
    lastEle = newEle;
 end
 
 %plot(y,concentration,'b*-')

end

function out = getSolution(variable,x,y,Ele,solution_ien,p_ien,vertexData,whichvariable)
%[ien_mesh,IENall,pAll] = elementIEN(Ele,meshData,solution_ien,p_ien);
[IENall,pAll,xy] = elementIEN(Ele,solution_ien,p_ien,vertexData);
%[JInverse, detJ,gij,xCord,yCord,A,hA] =  ...
%                    getjacobian(ien_mesh,vertexData,0,0);
     
 xCord = xy(1,:);
 yCord = xy(2,:);               
EvalPoints = getPoints(x,y,xCord,yCord);
%pAll(1)
%[ ShapeFunc, divShapeFunc ] = ShapeFunc_quad(EvalPoints,pAll(1));
p=pAll(1);
quadsf = QuadShapeFunc(EvalPoints,p);
[ShapeFunc, divShapeFunc] = quadsf.QuadShapeFuncTable();        
% solution = zeros(length(IENall),1);
% for i = 1:length(IENall)
%     solution(i) = u(IENall(i));
% end
nssl = length(IENall);
variable_ele = zeros(nssl,1);

for i = 1:nssl
    variable_ele(i) = variable(IENall(i));
end

%variable_ele
[u,v,p] = getuvp(variable_ele,pAll(1));
    
if strcmp(whichvariable ,'u')==1
        out = u'*ShapeFunc;
    elseif strcmp(whichvariable,'v')==1
       out = v'*ShapeFunc;
    elseif strcmp(whichvariable,'uv')==1
        out = sqrt((u'*ShapeFunc)^2+(v'*ShapeFunc)^2);
    elseif strcmp(whichvariable,'vorticity')==1
        out = -gradNx*v+gradNy*u;
    else
       out = p'*ShapeFunc;
end
end

function eleSearch = getSurrondEle(lastEle,meshData,solution_ien,p_ien,vertexData)
[ien_mesh] = elementIEN(lastEle,meshData,solution_ien,p_ien);
eleSearch = [];
for i = 1:3
    SurroundEle = vertexData{ien_mesh(i),3};
    eleSearch = [eleSearch,SurroundEle];
end
end

function [searchELE]=findNextEle(x,y,searchPatch,solution_ien,p_ien,vertexData)
%% search for the element contains x y point
distanceToCenter = zeros(length(searchPatch),1);
searchELE = 0;
for i = 1:length(searchPatch)
    ele = searchPatch(i);
    %[ien_mesh,IENall,pAll] = elementIEN(ele,meshData,solution_ien,p_ien);
    [IENall,pAll,xy] = elementIEN(ele,solution_ien,p_ien,vertexData);
     xCord = xy(1,:);
    yCord = xy(2,:);
    %[JInverse, detJ,gij,xCord,yCord,A,hA] =  ...
     %               getjacobian(ien_mesh,vertexData,0,0);
     inout = inpolygon(x,y,[xCord(1),xCord(2),xCord(4),xCord(3)],[yCord(1),yCord(2),yCord(4),yCord(3)]);
     if inout
         searchELE = ele;
     end
    %[EvalPoints, distanceToCenter(i)]=getPoints(x,y,xCord,yCord);

end
end


function [EvalPoints,distanceToCenter]=getPoints(x,y,xCordM1,yCordM1)
%% get center coordinate 
x1 = xCordM1(1);
x2 = xCordM1(2);
x3 = xCordM1(3);
x4 = xCordM1(4);

y1 = yCordM1(1);
y2 = yCordM1(2);
y3 = yCordM1(3);
y4 = yCordM1(4);

xmid = (x1+x2+x3+x4)/4;
ymid = (y1+y2+y3+y4)/4;

distanceToCenter = sqrt((x-xmid)^2+(y-ymid)^2);

%xi = (x-x1)/abs(x1-x2);
%eta = (y-y1)/abs(y4-y1);
%EvalPoints(1) = xi;
%EvalPoints(2) = eta;

lhs = [1,x1,y1,x1*y1;1,x2,y2,x2*y2;1,x3,y3,x3*y3;1,x4,y4,x4*y4];
rhs1 = [0;1;0;1];
rhs2 = [0;0;1;1];
coef1 = lhs\rhs1;
coef2 = lhs\rhs2;
EvalPoints(1) = [1,x,y,x*y]*coef1;
EvalPoints(2) = [1,x,y,x*y]*coef2;
% if xi>1
%     xi
%     x
%     y
%     xCordM1
%     yCordM1
% end

% %% get evalPoints in xi,eta coordinate
% a = (y4-y1)/(x4-x1);
% b = y1-a*x1;
% l1 = abs(a*x-y+b)/sqrt(a^2+1);
% l = sqrt((x1-x2)^2+(y2-y1)^2);
% xi = l1/l;
% a = (y2-y1)/(x2-x1);
% b = y1-a*x1;
% m1 = abs(a*x-y+b)/sqrt(a^2+1);
% m = sqrt((x1-x4)^2+(y4-y1)^2);
% eta = m1/m;
% EvalPoints(1) = xi;
% EvalPoints(2) = eta;
end