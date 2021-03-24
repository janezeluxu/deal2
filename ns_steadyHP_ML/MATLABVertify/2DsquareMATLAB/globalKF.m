function [lhs, rhs] = globalKF(iter,variable,totalDOF,meshData,IEN,vertexData,IBC)
global miu;
global Grid_size;
lhs = zeros(totalDOF,totalDOF);
rhs = zeros(totalDOF,1);
nInt = 3;
for ele = 1:Grid_size 
    ele;
    [IENall,xcord,ycord] = elementIEN(ele,meshData,vertexData,IEN);
    nssl = length(IENall);
    variable_ele = zeros(nssl,1);
    for i = 1:nssl
        variable_ele(i) = variable(IENall(i));
    end
    [k,f] = ElementMatrix(xcord,ycord,variable_ele,miu,nInt);
    %assemble    
    for i = 1:length(f)
        for j = 1:length(f)
            lhs(IENall(i),IENall(j)) = k(i,j)+lhs(IENall(i),IENall(j));
        end
        rhs(IENall(i)) = f(i)+rhs(IENall(i));
    end
end

[lhs,rhs] = applyBC(IBC,rhs,lhs);
end

function [lhs,rhs] = applyBC(IBC,rhs,lhs)
%% apply Strong boundary condition
for i = 1:length(IBC)
    rhs(IBC(i))=0;
    lhs(IBC(i),:) = 0;
    lhs(:,IBC(i)) = 0;
    lhs(IBC(i),IBC(i)) = 1;
end
% rhs(1) = 0;
% rhs(nPoints) = 0;
% 
% lhs(1,:) = 0;
% lhs(:,1) = 0;
% lhs(1,1) = 1;
% lhs(nPoints,:) = 0;
% lhs(:,nPoints) = 0;
% lhs(nPoints,nPoints) = 1;
end
