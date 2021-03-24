function [] = writevtkfile(meshData,vertexData,u,p)
fileID = fopen('./testvtk/elements.txt','w');
formatSpec = '%d %d %d %d\n';
[nrows,ncols] = size(meshData);
for row = 1:nrows
    fprintf(fileID,formatSpec,meshData(row,:));
end
fclose(fileID);

fileID = fopen('./testvtk/nodes.txt','w');
formatSpec = '%2.16f %2.16f \n';
[nrows,ncols] = size(vertexData);
for row = 1:nrows
    fprintf(fileID,formatSpec,vertexData(row,:));
end

fclose(fileID);

fileID = fopen('./testvtk/pointdata.txt','w');
formatSpec = '%2.16f %2.16f %2.16f \n';
[nrows,ncols] = size(u);
value = zeros(nrows,3);
value(:,1) = u(:,1);
value(:,2) = u(:,2);
value(:,3) = p;
for row = 1:nrows
    fprintf(fileID,formatSpec,value(row,:));
end
fclose(fileID);

% fileID = fopen('./testvtk/cases/solution/nonUniform/pointdataStatic.txt','w');
% formatSpec = '%2.16f %2.16f %2.16f \n';
% [nrows,ncols] = size(u);
% value = zeros(nrows,3);
% value(:,1) = u(:,1);
% value(:,2) = u(:,2);
% value(:,3) = p;
% for row = 1:nrows
%     fprintf(fileID,formatSpec,value(row,:));
% end
% fclose(fileID);

% fileID = fopen('./testvtk/celldata.txt','w');
% formatSpec = '%2.16f \n';
% [nrows,ncols] = size(tau);
% value = zeros(nrows,1);
% value(:,1) = tau;
% %value(:,2) = tau;
% %value(:,3) = gGlobal;
% for row = 1:nrows
%     fprintf(fileID,formatSpec,value(row,:));
% end
% fclose(fileID);
end