function [IBCuv,IBCp,BC_val_uv,BC_val_p] = readin_bc(bcfileName,totalNode)
%bcfileName = "../mesh/p2bc.txt";
fileID = fopen(bcfileName,'r');
tline = fgetl(fileID);
i = 1;
while ischar(tline)    
    Array{i}=tline;
    i=i+1;
    tline = fgetl(fileID);
end
fclose(fileID);

IBC = zeros(length(Array),1);
BC_val = zeros(length(Array),1);
for i = 1:length(Array)
    array1 = split(Array{i},'=');
    IBC(i) = str2num(array1{1})+1;
    BC_val(i) = str2num(array1{2});
end


[Indexp] = find(IBC>totalNode);
IBCp = IBC(Indexp)-totalNode;
BC_val_p = BC_val(Indexp);

[Indexuv] = find(IBC<totalNode);
IBCuv = IBC(Indexuv);
BC_val_uv = BC_val(Indexuv);

% [Indexu] = find(IBC<totalNode & mod(IBC,2));
% IBCu = (IBC(Indexu)+1)/2;
% BC_val_u = BC_val(Indexu);
% 
% [Indexv] = find(IBC<totalNode & ~mod(IBC,2));
% IBCv = IBC(Indexv)/2;
% BC_val_v = BC_val(Indexv);
end