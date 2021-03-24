function []=getu2()
clc
clear all

mesh_file1 = "../mesh/p1mesh.txt";
mesh_file2 = "../mesh/p3mesh.txt";
bcfileName = "../mesh/p3bc.txt";

uv = load("../solution/100-un-3191.txt");
p = load("../solution/100-pn-3191.txt");
p_t = load("../solution/100-pt_n-3191.txt");
uv_t = load("../solution/100-ut_n-3191.txt");

totalNode1 = length(uv);
%% setup mesh
p1 = 1;
p2 = 3; 
disp('load mesh_file1')
[IEN_all1,p_All1,vertexData1] = readin_mesh(mesh_file1);
Grid_size = length(p_All1);
IEN_uv1 = zeros(Grid_size,2*(p1+1)^2);
IEN_p1 = zeros(Grid_size,(p1+1)^2);
for ele = 1:Grid_size
    [uIEN,vIEN,uvIEN,pIEN] = getuvp(IEN_all1(ele,:),p_All1(ele));
    IEN = pIEN-totalNode1;
    IEN_p1(ele,:) = IEN;
    IEN_uv1(ele,:) = uvIEN;
end
   
disp('load mesh_file2')
[IEN_all2,p_All2,vertexData2] = readin_mesh(mesh_file2);
%IEN_all2(1,:)
Grid_size2 = length(p_All2);
IEN_uv2 = zeros(Grid_size2,2*(p2+1)^2);
IEN_p2 = zeros(Grid_size2,(p2+1)^2);
totalNode2 = IEN_all2(1,3)-1;
for ele = 1:Grid_size2
    [uIEN,vIEN,uvIEN,pIEN] = getuvp(IEN_all2(ele,:),p_All2(ele));
    IEN = pIEN-totalNode2;
    IEN_p2(ele,:) = IEN;
    IEN_uv2(ele,:) = uvIEN;
end
%IEN_uv2(1,:)

[IBCuv2,IBCp2,BC_val_uv2,BC_val_p2] = readin_bc(bcfileName,totalNode2);

disp('solve uv2')
[velocity2]=mapuv(IBCuv2,BC_val_uv2,p2,IEN_uv2,p_All2,vertexData2,p1,IEN_uv1,p_All1,vertexData1,uv);
%size(u2)

disp('solve uvt2')
[velocity2_t]=mapuv(IBCuv2,BC_val_uv2,p2,IEN_uv2,p_All2,vertexData2,p1,IEN_uv1,p_All1,vertexData1,uv_t);

 disp('solve p2')
 [pressure2]=mapSolution(IBCp2,BC_val_p2,p2,IEN_p2,p_All2,vertexData2,p1,IEN_p1,p_All1,vertexData1,p);

 disp('solve pt2')
 [pressure2_t]=mapSolution(IBCp2,BC_val_p2,p2,IEN_p2,p_All2,vertexData2,p1,IEN_p1,p_All1,vertexData1,p_t);

%% write to file
filename = strcat('../solution/100-un-3192.txt');
fileID = fopen(filename,'w');
fprintf(fileID,'%6.20f ',velocity2);
fclose(fileID);

filename = strcat("../solution/100-ut_n-3192.txt");
fileID = fopen(filename,'w');
fprintf(fileID,'%6.20f ',velocity2_t);
fclose(fileID);

filename = strcat("../solution/100-pn-3192.txt");
fileID = fopen(filename,'w');
fprintf(fileID,'%6.20f ',pressure2);
fclose(fileID);

filename = strcat("../solution/100-pt_n-3192.txt");
fileID = fopen(filename,'w');
fprintf(fileID,'%6.20f ',pressure2_t);
fclose(fileID);

plotSolution(velocity2',pressure2',mesh_file2)
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