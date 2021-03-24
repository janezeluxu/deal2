function [eleBC] = readin_bcele(bcfileName)
%bcfileName = "../mesh/bcEle.txt";
fileID = fopen(bcfileName,'r');
tline = fgetl(fileID);
i = 1;
while ischar(tline)    
    Array{i}=tline;
    i=i+1;
    tline = fgetl(fileID);
end
fclose(fileID);

eleBC = zeros(length(Array),3);
for i = 1:length(Array)
    eleBC(i,:) = str2num(Array{i});
end
end