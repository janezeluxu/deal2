function [F] = calculateForce(variable,BCF,lam_BCF,solution_ien,p_ien,order,vertexData)
%global nsd;
%global miu;
F = [0;0];
miu = 0.01;
%ndirec = {@(x,y) (2*(x-10));@(x,y) (2*(y-10))};
%nFdirec = [1 0 -1 0; 0 -1 0 1];
ep = 1e-3;
%% loop over element
for i = 1:length(BCF)
%for i = 1:1 
    ele = BCF(i);
    lam = lam_BCF(i);
    
    [IENall,pAll,xy] = elementIEN(ele,solution_ien,p_ien,vertexData);
    nssl = length(IENall);
    variable_ele = zeros(nssl,1);
    
    for j = 1:nssl
        variable_ele(j) = variable(IENall(j));
    end
    [u,v,uv,p] = getuvp(variable_ele,pAll);
    
    xCord = xy(1,:);
    yCord = xy(2,:);
    
    n = nIntergerPoints(max(pAll),0);
    [xi,w] =  GaussQuad(n, 1);
    
    [h,dirxy] = getBCedge(xCord,yCord,lam);
    qPoints = getBCQuadrature(xi,w,lam);
    
    quadsf = QuadShapeFunc(qPoints,order);
    [ShapeFunc, divShapeFunc] = quadsf.QuadShapeFuncTable();
     
    quadsf1 = QuadShapeFunc(qPoints,1);
    [ShapeFunc1, ~] = quadsf1.QuadShapeFuncTable();
     x_dir = ones(1,4)*dirxy(1);
     y_dir = ones(1,4)*dirxy(2);
%     ShapeFunc1
    Jw = h*qPoints(:,3);
    nP = size(qPoints,1);
    %% loop over quadrature points
    for k = 1:nP 
        %k
        xi = qPoints(k,1);
        eta = qPoints(k,2);
        [JInverse, detJ,gij] = getjacobian(xCord,yCord,xi,eta);
        
        gradNGlobal = (divShapeFunc(:,:,k)*JInverse)';
        gradNx = gradNGlobal(1,:);
        gradNy = gradNGlobal(2,:);
                
        
        ux = gradNx*u;
        uy = gradNy*u;
         
        vx = gradNx*v;
        vy = gradNy*v;
        
        pquad = (ShapeFunc(:,k))'*p;
        
        ndir = (ShapeFunc1(:,k))'*[x_dir;y_dir]';
        
        shereStress = [miu*(ux+ux),miu*(ux+vy);miu*(ux+vy),miu*(vy+vy)];
        
        F = F-shereStress*ndir'*Jw(k)-(pquad*ndir')*Jw(k);
        %F = F+shereStress*ndir'*Jw(k);
        %L = L+(pquad*yydir)*Jw(k);       
    end
end

end

function [u,v,uv,p] = getuvp(variable_ele,pAll)
uvsize = 3*(pAll+1)^2;
if pAll <3
    u = variable_ele(1:3:uvsize);
    v = variable_ele(2:3:uvsize);
    p = variable_ele(3:3:end);
    uv = variable_ele([1:3:uvsize,2:3:uvsize]);
elseif pAll==3
    %1:3:uvsize
    u = variable_ele([1,4,7,10,13,14,19,20,25,26,31,32,37,38,39,40]);
    v = variable_ele([2,5,8,11,15,16,21,22,27,28,33,34,41,42,43,44]);
    p = variable_ele([3,6,9,12,17,18,23,24,29,30,35,36,45,46,47,48]);
    uv = variable_ele([1,4,7,10,13,14,19,20,25,26,31,32,37,38,39,40,...
        2,5,8,11,15,16,21,22,27,28,33,34,41,42,43,44]);
elseif pAll==5
    %1:3:uvsize
    u = variable_ele([1,4,7,10,13:16,25:28,37:40,49:52,61:76]);
    v = variable_ele([2,5,8,11,17:20,29:32,41:44,53:56,77:92]);
    p = variable_ele([3,6,9,12,21:24,33:36,45:48,57:60,93:108]);
    uv = variable_ele([1,4,7,10,13:16,25:28,37:40,49:52,61:76,...
        2,5,8,11,17:20,29:32,41:44,53:56,77:92]);
 elseif pAll==7
    u = variable_ele([1,4,7,10,13:18,31:36,49:54,67:72,85:120]);
    v = variable_ele([2,5,8,11,19:24,37:42,55:60,73:78,121:156]);
    p = variable_ele([3,6,9,12,25:30,43:48,61:66,79:84,157:192]);
    uv = variable_ele([1,4,7,10,13:18,31:36,49:54,67:72,85:120,...
        2,5,8,11,19:24,37:42,55:60,73:78,121:156]);   
end
end

function [h,dir] = getBCedge(xCord,yCord,edgeNumber)
if (edgeNumber ==1)
    % distance between v1 and v3
    h = ((xCord(1)-xCord(3))^2+(yCord(1)-yCord(3))^2)^0.5;
    xx = xCord(1)-xCord(3);
    yy = yCord(1)-yCord(3);
    
    x1 = xCord(1);
    y1 = yCord(1);
elseif (edgeNumber ==2)
    % distance between v2 and v4
    h = ((xCord(2)-xCord(4))^2+(yCord(2)-yCord(4))^2)^0.5;
    xx = xCord(2)-xCord(4);
    yy = yCord(2)-yCord(4);
    
    x1 = xCord(2);
    y1 = yCord(2);
elseif edgeNumber ==3
    % distance between v3 and v4
    h = ((xCord(3)-xCord(4))^2+(yCord(3)-yCord(4))^2)^0.5;
    xx = xCord(3)-xCord(4);
    yy = yCord(3)-yCord(4);
    
    x1 = xCord(3);
    y1 = yCord(3);
elseif (edgeNumber ==4)
    % distance between v1 and v2
    h = ((xCord(1)-xCord(2))^2+(yCord(1)-yCord(2))^2)^0.5;
    xx = xCord(1)-xCord(4);
    yy = yCord(1)-yCord(4);
    
    x1 = xCord(1);
    y1 = yCord(1);
end

    if xx ==0
        if (x1-9.5) ==0
            dir = [-1,0];
        else
            dir = [1,0];
        end
    end
    if yy ==0
        if (y1-9.5) ==0
            dir = [0,-1];
        else
            dir = [0,1];
        end
    end
end

function qPoints = getBCQuadrature(xi,w,edgeNumber)
if (edgeNumber ==1)
    lxi = length(xi);
    qPoints = zeros(lxi,3);
    for i = 1:lxi
        qPoints(i,:)=[0,xi(i),w(i)];
    end
elseif (edgeNumber ==2)
    lxi = length(xi);
    qPoints = zeros(lxi,3);
    for i = 1:lxi
        qPoints(i,:)=[1,xi(i),w(i)];
    end
elseif (edgeNumber ==3)
    lxi = length(xi);
    qPoints = zeros(lxi,3);
    for i = 1:lxi
        qPoints(i,:)=[xi(i),1,w(i)];
    end
elseif (edgeNumber ==4)
    lxi = length(xi);
    qPoints = zeros(lxi,3);
    for i = 1:lxi
        qPoints(i,:)=[xi(i),0,w(i)];
    end
end
end
    