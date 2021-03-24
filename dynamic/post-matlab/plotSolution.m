function plotSolution(uv,p,mesh_file)
%clc
%clear all
%uv = load("../solution/100-un-3192.txt");
%p = load("../solution/100-pn-3192.txt");
%mesh_file = "../mesh/p2mesh.txt";
variable = [uv,p];
% length(uv)
% variable(1)
% variable(2)
% variable(191961)
% p(1)
[solution_ien,p_ien,vertexData] = readin_mesh(mesh_file);
Grid_size = size(solution_ien,1);
polyplot(Grid_size,solution_ien,p_ien,vertexData,variable);
end