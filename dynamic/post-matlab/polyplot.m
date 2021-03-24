function polyplot(Grid_size,solution_ien,p_ien,vertexData,variable)
%global TotalDOF; 
%create uniform p shape function tables, p is a list
for ele = 1:Grid_size     
%for ele = 1:1
    %ele
% element = [1,7];
% for in = 1:length(element)
%     ele = element(in);
    %if uniform p, use tabulated shapefunction and div values
    %else, calculate mixorder element shape functions  
    %% get meshIEN, solutionIEN, and shape functions
    [IENall,pAll,xy] = elementIEN(ele,solution_ien,p_ien,vertexData);
    nssl = length(IENall);
    %edgeU = edge_usage(ele,:);
    variable_ele = zeros(nssl,1);
    
    for i = 1:nssl
        variable_ele(i) = variable(IENall(i));
    end
    
    %variable_ele'
    %variable_ele
    [u,v,p] = getuvp(variable_ele,pAll);
    %IENall
    x = xy(1,:);
    y = xy(2,:);
    %% plot high order solution for each element 
    dx = 0:0.2:1;
    [xc,yc,uIntPoint,sx,sy,uexact] = ElementSurfacePlot(pAll,u,dx,x,y);
    
    %subplot(2,2,plotnum)
    sp = surf(xc,yc,uIntPoint);
    %set(sp,'EdgeColor',[1 1 1]);
    %set(sp,'FaceColor',[1 1 1]);
    %colormap default
    sp.EdgeColor = 'none';
    
    %sp = contour(xc,yc,uIntPoint,[0.1:0.2:1.0],'linewidth',2);
    %sp = contour(xc,yc,uexact,[0.1:0.2:1.0],'linewidth',2);
    %uIntPoint
    %sp = contour(sx,sy,uIntPoint,[0.1:0.2:1.0],'linewidth',2);
    %sp = contour(xc,yc,uIntPoint,[-0.2:0.04:1.2]);
    hold on
    
end
%light
%light('Position',[0.75 0.5 1.2],'Style','infinite')
end

function [u,v,p] = getuvp(variable_ele,pAll)
%pAll
uvsize = 3*(pAll+1)^2;
if pAll <3
    u = variable_ele(1:3:uvsize);
    v = variable_ele(2:3:uvsize);
    p = variable_ele(3:3:end);
elseif pAll==3
    %1:3:uvsize
    u = variable_ele([1,4,7,10,13,14,19,20,25,26,31,32,37,38,39,40]);
    v = variable_ele([2,5,8,11,15,16,21,22,27,28,33,34,41,42,43,44]);
    p = variable_ele([3,6,9,12,17,18,23,24,29,30,35,36,45,46,47,48]);   
elseif pAll==5
    %1:3:uvsize
    u = variable_ele([1,4,7,10,13:16,25:28,37:40,49:52,61:76]);
    v = variable_ele([2,5,8,11,17:20,29:32,41:44,53:56,77:92]);
    p = variable_ele([3,6,9,12,21:24,33:36,45:48,57:60,93:108]);
elseif pAll==7
    u = variable_ele([1,4,7,10,13:18,31:36,49:54,67:72,85:120]);
    v = variable_ele([2,5,8,11,19:24,37:42,55:60,73:78,121:156]);
    p = variable_ele([3,6,9,12,25:30,43:48,61:66,79:84,157:192]);
    uv = variable_ele([1,4,7,10,13:18,31:36,49:54,67:72,85:120,...
        2,5,8,11,19:24,37:42,55:60,73:78,121:156]);   
end
end

function [xc,yc,uIntPoint,sx,sy,uexact] = ElementSurfacePlot(pAll,variable_ele,dx,x,y)

[xi,eta] = meshgrid(dx,dx);

uIntPoint = zeros(size(xi,1),size(eta,2));
uexact = zeros(size(xi,1),size(eta,2));
xc = zeros(size(xi,1),size(eta,2));
yc = zeros(size(xi,1),size(eta,2));
sx = zeros(size(xi,1),size(eta,2));
sy = zeros(size(xi,1),size(eta,2));
for i = 1:size(xi,1)
    for j = 1:size(eta,1)
        
        N1 = (1-xi(i,j))*(1-eta(i,j));
        N2 = xi(i,j)*(1-eta(i,j));
        N3 = xi(i,j)*eta(i,j);
        N4 = (1-xi(i,j))*eta(i,j);
        
        xc(i,j) = N1*x(1)+N2*x(2)+N3*x(4)+N4*x(3);
        yc(i,j) = N1*y(1)+N2*y(2)+N3*y(4)+N4*y(3);
        
        sx(i,j) = (xc(i,j)+yc(i,j))/(sqrt(2));
        sy(i,j) = abs(xc(i,j)-yc(i,j))/(sqrt(2));
        
        qPoints(1) = xi(i,j);
        qPoints(2) = eta(i,j);
        
        if range(pAll) == 0
            %disp('uniform p element');
            p=pAll(1);
            quadsf = QuadShapeFunc(qPoints,p);
            [ShapeFunc, divShapeFunc] = quadsf.QuadShapeFuncTable();
            
        else
            %disp('nonuniform p element');
            %simplexsf = SimplexShapeFunc(qPoints,IENall,pAll);
            %[ShapeFunc, ~] = simplexsf.variablePShapeFunctionTable();
            
        end
        
        uIntPoint(i,j) = variable_ele'*ShapeFunc;
        
        %uexact(i,j) = Exact(xc(i,j),yc(i,j),casenumber,a,kappa);
    end
end
end