function [Grid_size,vertices_size,Edge_size,total_DOF,mesh_IEN,...
    mesh_Edge,Edge_v,IBC,IBC_edge] ...
    = load_mesh(fileName,bcFile)
%fileName = "../output/mesh/Shape1-1EdgeInfo.txt";
%bcFile = "../output/mesh/Shape1-1bc.txt";
fileID = fopen(fileName,'r');
tline = fgetl(fileID);
i = 1;
cellcount = 1;
while ischar(tline)
    C = textscan(tline,'%s');
    c= C{1};
    if ( strcmp(c{1},'Grid_size')==1)
         Grid_size = str2num(c{2});
    end
    if ( strcmp(c{1},'vertices_size')==1)
         vertices_size = str2num(c{2});
    end
    if ( strcmp(c{1},'Edge_size')==1)
         Edge_size = str2num(c{2});
    end
    if ( strcmp(c{1},'total_DOF')==1)
         total_DOF = str2num(c{2});
    end
    if ( strcmp(c{1},'Face_Vertex_Edge')==1)
         startline = i;
    end
    if ( strcmp(c{1},'Edge_Vertex')==1)
         startlineEdge = i;
    end
    Array{i}=tline;
    i=i+1;
    tline = fgetl(fileID);
end
fclose(fileID);

mesh_IEN = zeros(Grid_size,4);
mesh_Edge = zeros(Grid_size,4);
for n = 1:Grid_size
    mesh_IEN(n,:) = str2num(Array{startline+3*(n-1)+2})+1;
    mesh_Edge(n,:) = str2num(Array{startline+3*(n-1)+3})+1;
end

Edge_v = zeros(Edge_size,2);
%for n = 1:Edge_size
%    Edge_v(n,:) = str2num(Array{startlineEdge+n})+1;
%end

fileID = fopen(bcFile,'r');
tline = fgetl(fileID);
i = 1;
count = 0;
while ischar(tline)
    C = textscan(tline,'%s');
    %c= C{1};
    Array{i}=tline;
    i=i+1;
    tline = fgetl(fileID);
    count = count+1;
    %BCb(count,:) = str2num(Array{startlineEdge+n})+1;
end
fclose(fileID);

BCb = zeros(count,6);
for i = 1:count
    BCb(i,:) = str2num(Array{i})+1;
end

IBC_vertex = zeros(2*count,2);
for i = 1:count
    IBC_vertex(2*i-1,1) = BCb(i,5);
    IBC_vertex(2*i,1) = BCb(i,6);
    IBC_vertex(2*i-1,2) = BCb(i,2);
    IBC_vertex(2*i,2) = BCb(i,2);
end
[C,idx]=unique(IBC_vertex(:,1),'stable');
IBC = IBC_vertex(idx,:);

IBC_edge = zeros(count,1);
for i = 1:count
    IBC_edge(i,1) = BCb(i,4);
    IBC_edge(i,2) = BCb(i,2);
end

k = IBC(:,1) == 1;
IBC(k,2) = 5;

k = IBC(:,1) == 2;
IBC(k,2) = 6;

k = IBC(:,1) == 3;
IBC(k,2) = 7;

k = IBC(:,1) == 4;
IBC(k,2) = 8;
end