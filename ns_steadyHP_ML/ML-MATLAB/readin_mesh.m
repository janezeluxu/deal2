function [solution_ien,p_ien,cellvertexDatax,cellvertexDatay,vertexData,v_e] ...
    = readin_mesh(fileName,fileNameVertex)
%fileName = "../output/mesh/lidcavity-1mesh.txt";
%fileName
fileID = fopen(fileName,'r');
tline = fgetl(fileID);
i = 1;
cellcount = 1;
while ischar(tline)
    C = textscan(tline,'%s');
    c= C{1};
    if ( strcmp(c{1},'cellnum')==1)
         startline(cellcount) = i;
         cellnum(cellcount) = str2num(c{2});
         nIEN(cellcount) = str2num(c{3});
         cellcount = cellcount+1;
    end
    
    Array{i}=tline;
    i=i+1;
    tline = fgetl(fileID);
end
fclose(fileID);

for n = 1:length(startline)
    nLine = startline(n);
    for m = 1:nIEN(n)
        solution_ien(n,m) = str2num(Array{nLine+m})+1;
    end
    p_ien(n) = sqrt(nIEN(n)/3)-1;
    cellvertexDatax(n,:) = str2num(Array{nLine+m+1});
    cellvertexDatay(n,:) = str2num(Array{nLine+m+2});
end

%cellvertexData{1} = cellvertexDatax;
%cellvertexData{2} = cellvertexDatay;
%cellvertexData

%solution_ien
%vertexDatax
%vertexDatay

%fileNameVertex = "../output/mesh/lidcavity-1v_to_e_indices.txt";

fileID = fopen(fileNameVertex,'r');
tline = fgetl(fileID);
i = 1;
v_size = 0;
while ischar(tline)
    C = textscan(tline,'%s');
    c= C{1};
    Array{i}=tline;
    i=i+1;
    tline = fgetl(fileID);
    v_size = v_size+1;
    %tline
end
fclose(fileID);

vertex_size = v_size-1;
%grid_size = size(solution_ien,1);
%vertex_size = 625;
vertexData = zeros(vertex_size,2);
v_e = cell(vertex_size,1);
for v = 1:vertex_size
    list = str2num(Array{v});
    vertexData(v,:) = list(2:3);
    v_e{v,1} = list(4:end)+1; % start from 1
    %vertexData(v,:) = str2num(Array{v});
end
%solution_ien
%vertexData
%v_e
end