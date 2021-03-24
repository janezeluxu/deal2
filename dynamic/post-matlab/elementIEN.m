function [IENall,pAll,xy] = elementIEN(ele,solution_ien,p_ien,vertexData)
%get the IEN array and P list for a given element
%IEN_mesh = ien_map(ele,:);
IENall = solution_ien(ele,:);
IENall = IENall(IENall~=0);
pAll = p_ien(ele);
x = vertexData{1}(ele,:);
y = vertexData{2}(ele,:);
xy(1,:) = x;
xy(2,:) = y;
end