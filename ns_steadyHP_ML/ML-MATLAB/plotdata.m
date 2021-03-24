function plotdata()
global Grid_size
fileName = "../output/mesh/lidcavity_vms_result/lidcavity-800-5mesh.txt";
fileNameVertex = "../output/mesh/lidcavity_vms_result/lidcavity-800-5v_to_e_indices.txt";
%cellLevel = load("prediction_Shape1_Re50");
IENfilename = "../output/mesh/lidcavity_vms_result/lidcavity-800-5EdgeInfo.txt";
bcFile = "../output/mesh/lidcavity_vms_result/lidcavity-800-5bc.txt";

%outmeshfile = "../NS-steadyHO-quad/mesh/Shape1.txt";
%outbcfile = "../NS-steadyHO-quad/mesh/Shape1bc.txt";
%solutionFile = "../NS-steadyHO-quad/testvtk/cases/solution/cylinder/shape1.txt";
[solution_ien,p_ien,cellvertexDatax,cellvertexDatay,vertexData,v_e] ...
    = readin_mesh(fileName,fileNameVertex);
[Grid_size,n] = size(solution_ien)

[Grid_size,vertices_size,Edge_size,total_DOF,meshData,mesh_Edge,Edge_v,IBC,IBC_edge] ...
    = load_mesh(IENfilename,bcFile);

uv = load("../output/solution/lidcavity_vms_result/lidcavityRe800un-5-.txt");
p = load("../output/solution/lidcavity_vms_result/lidcavityRe800pn-5-.txt");

solution = [uv,p];

dx = 0.1:0.1:0.9;
dy = 0.1:0.1:0.9;
[X,Y]=meshgrid(dx,dy);
sizeN = length(dx);
p_sol = zeros(sizeN,sizeN);
u_sol = zeros(sizeN,sizeN);
v_sol = zeros(sizeN,sizeN);

for i = 1:length(dx)
    y = Y(i,:);
    x = X(i,:);    
    [u,v,p]=LinePlot(meshData,vertexData,cellvertexDatax,cellvertexDatay,...
        solution_ien,p_ien,solution,x,y);
    u_sol(i,:) = u;
    v_sol(i,:) = v;
    p_sol(i,:) = p;
    %u(i,:) = LinePlot(meshData,vertexData,solution_ien,p_ien,solution,x,y);
    %v(i,:) = LinePlot(meshData,vertexData,solution_ien,p_ien,solution,x,y);
end

%var

%var
figure(1)
contour(X,Y,u_sol,'linewidth',2)
figure(2)
contour(X,Y,v_sol,'linewidth',2)
figure(3)
contour(X,Y,p_sol,'linewidth',2)
end

function [u_sol,v_sol,p_sol]=LinePlot(meshData,vertexData,cellvertexDatax,cellvertexDatay,solution_ien,...
    order_list,var,x,y)
 global Grid_size;

 %u = load(solution);
 
 u_sol = zeros(1,length(y));
 v_sol = zeros(1,length(y));
 p_sol = zeros(1,length(y));
 x(1)
 y(1)
 %% start the search from boundary elements
 [startingEle]=findNextEle(x(1),y(1),1:Grid_size,...
        cellvertexDatax,cellvertexDatay);
 %startingEle = 12461;
 lastEle = startingEle;
 
 for i = 1:length(x)
     %lastEle
     xx = x(i)
     yy = y(i)
    %eleSearch = getSurrondEle(lastEle,meshData,...
    %    solution_ien,order_list,vertexData);
    %[newEle] = findNextEle(x(i),y(i),eleSearch,...
    %    cellvertexDatax,cellvertexDatay);
    
    %if newEle ==0
        [newEle] = findNextEle(x(i),y(i),1:Grid_size,...
        cellvertexDatax,cellvertexDatay);
    %end
    %x(i),y(i)
    newEle
    [u,v,p] = getSolution(var,x(i),y(i),newEle,meshData,...
        solution_ien,order_list,cellvertexDatax,cellvertexDatay);
    %cc = concentration(i)
    u_sol(i) = u
    v_sol(i) = v
    p_sol(i) = p
    %lastEle = newEle;
 end
 
 %plot(y,concentration,'b*-')
end

function [u,v,p] = getSolution(solution,x,y,ele,meshData,solution_ien,p_ien,cellvertexDatax,cellvertexDatay)
%[ien_mesh,IENall,pAll] = elementIEN(ele,meshData,solution_ien,p_ien);
%ien_mesh = meshData(ele,:);
IENall = solution_ien(ele,:);
IENall = IENall(IENall~=0);
pAll = p_ien(ele);

xCord = cellvertexDatax(ele,:);
yCord = cellvertexDatay(ele,:);
        
EvalPoints = getPoints(x,y,xCord,yCord);
%pAll(1)
[ ShapeFunc, divShapeFunc ] = ShapeFunc_quad(EvalPoints,pAll(1));

%solutionsize = size(solution)
variable_ele = zeros(length(IENall),1);
for i = 1:length(IENall)
    variable_ele(i) = solution(IENall(i));
end
%ShapeFunc
%variable_ele
nssl = length(IENall);
nsd = 2;
nsf = nssl/(nsd+1); 
    
u_ele = variable_ele(1:nsf);
v_ele = variable_ele(nsf+1:2*nsf);
p_ele = variable_ele(2*nsf+1:3*nsf);

u = u_ele'*ShapeFunc;
v = v_ele'*ShapeFunc;
p = p_ele'*ShapeFunc;
end

function eleSearch = getSurrondEle(lastEle,meshData,solution_ien,p_ien,vertexData)
%[ien_mesh] = elementIEN(lastEle,meshData,solution_ien,p_ien);
ien_mesh = meshData(lastEle,:);
%IENall = IENall(IENall~=0);
%pAll = p_ien(ele);

%xCord = cellvertexDatax(ele,:);
%yCord = cellvertexDatay(ele,:);

eleSearch = [];
for i = 1:3
    SurroundEle = vertexData{ien_mesh(i),3};
    eleSearch = [eleSearch,SurroundEle];
end
end

function [searchELE]=findNextEle(x,y,searchPatch,cellvertexDatax,cellvertexDatay)
%% search for the element contains x y point
%x
%y
%cellvertexDatax(1:100)
distanceToCenter = zeros(length(searchPatch),1);
searchELE = 0;
for i = 1:length(searchPatch)
    ele = searchPatch(i);
    %[ien_mesh,IENall,pAll] = elementIEN(ele,meshData,solution_ien,p_ien);
    %ien_mesh = meshData(ele,:);
    xCord = cellvertexDatax(ele,:);
    yCord = cellvertexDatay(ele,:);
    %[JInverse, detJ,gij,xCord,yCord,A,hA] =  ...
    %                getjacobian(ien_mesh,vertexData,0,0);
                
    in = inpolygon(x,y,[xCord(1),xCord(2),xCord(4),xCord(3)],[yCord(1),yCord(2),yCord(4),yCord(3)]);
    if in
       searchELE = ele;
    end
    
%     if i == 929
%         x,y,xCord,yCord
%         in
%     end
    
%     [EvalPoints, distanceToCenter(i)]=getPoints(x,y,xCord,yCord);
%     if EvalPoints(1)>=0 && EvalPoints(1)<=1 && EvalPoints(2)>=0 && EvalPoints(2)<=1
% %         x
% %         y
% %         xCord
% %         yCord
% %         ien_mesh
% %         EvalPoints
%         searchELE = ele
%         %searchELE
%     end
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

xi = (x-x1)/abs(x3-x1);
eta = (y-y2)/abs(y2-y1);
EvalPoints(1) = xi;
EvalPoints(2) = eta;

% if xi>1
%     xi
%     x
%     y
%     xCordM1
%     yCordM1
% end

end