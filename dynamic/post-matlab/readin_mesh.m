function [solution_ien,p_ien,vertexData] = readin_mesh(fileName)
%fileName = "../../mesh/1mesh.txt";
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
    vertexDatax(n,:) = str2num(Array{nLine+m+1});
    vertexDatay(n,:) = str2num(Array{nLine+m+2});
end

vertexData{1} = vertexDatax;
vertexData{2} = vertexDatay;
end
