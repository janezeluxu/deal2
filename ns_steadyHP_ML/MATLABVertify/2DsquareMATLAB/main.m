function main()
clc
clear all
global Newton_maxIter;
global miu;
global Grid_size;
global TotalNode;
miu = 1;
Newton_maxIter = 7;
N = 11;
M = 11;
[meshData,IEN,vertexData,IBC,BCval,TotalNode,totalDOF] = createMesh(N,M);
variable = starting(IBC,BCval,TotalNode);
[variable] = Newton_solve(meshData,IEN,vertexData,IBC,totalDOF,variable)

u = variable(1:TotalNode);
v = variable(TotalNode+1:2*TotalNode);
p = variable(2*TotalNode+1:3*TotalNode);

size(variable)
filename = 'solution.txt';
writeSolution(filename,variable);

refSolutionfile = 'refSolution.txt';
[L2_ref,H1_ref] = getErrorwithRef(variable,refSolutionfile,N,M);
[L2_est,H1_est] = estimateError(variable,meshData,IEN,vertexData);
L2_ref
L2_est
% eL2 = L2_ele./L2_ref;
% eH1 = H1_ele./H1_ref;
% 
norm(L2_est)
norm(L2_ref)
 global_L2e = norm(L2_est)/norm(L2_ref)
global_H1e = norm(H1_est)/norm(H1_ref)

%surfaceplot(u,v,p,M,N,eL2)
end

function [variable] = starting(IBC,BCval,TotalNode)
u = ones(TotalNode,1)*0;
v = zeros(TotalNode,1);
p = zeros(TotalNode,1);

variable = [u;v;p];

for i = 1:length(IBC)
    variable(IBC(i)) = BCval(i);
end

end

function writeSolution(filename,variable)
fileID = fopen(filename,'w');
formatSpec = '%2.16f \n';
fprintf(fileID,formatSpec,variable);
fclose(fileID);
end

function surfaceplot(u,v,p,M,N,eff)
figure(1)
[X,Y] = meshgrid(1:N, 1:M);
Z = zeros(M,N);
for i = 1:N
    Z(:,i) = u(1+(i-1)*M:i*M);
    
end
surf(X,Y,Z)
figure(2)
[X,Y] = meshgrid(1:N, 1:M);
Z = zeros(M,N);
for i = 1:N
    Z(:,i) = v(1+(i-1)*M:i*M);
    
end
surf(X,Y,Z)
figure(3)
[X,Y] = meshgrid(1:N, 1:M);
Z = zeros(M,N);
for i = 1:N
    Z(:,i) = p(1+(i-1)*M:i*M);
    
end
surf(X,Y,Z)

figure(4)
[X,Y] = meshgrid(1:N-1, 1:M-1);
Z = zeros(M-1,N-1);
for i = 1:N-1
    Z(:,i) = eff(1+(i-1)*(M-1):i*(M-1));
    
end
surf(X,Y,Z)
end