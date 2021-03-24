xcoord = load('x.mat');
x = xcoord.xx;
ycoord = load('y.mat');
y = ycoord.yy;
for i = 1:length(x)
    
end

outputfilename = 'shape2.txt'
fileID = fopen(outputfilename,'w');
%formatSpec = '%s %f \n';
formatSpec = 'Point(%d)={ %2.6f, %2.6f, 0, lc2};\n';
%fprintf(formatSpec,A1,A2)
%[nrows,ncols] = size(cellSize);
%value = zeros(nrows,1);
%value(:,1) = cellSize;
for i = 1:length(x)
    fprintf(fileID,formatSpec,100+i,x(i),y(i));
    %formatSpec = 'X is %4.2f meters or %8.3f mm\n';
end
fclose(fileID);