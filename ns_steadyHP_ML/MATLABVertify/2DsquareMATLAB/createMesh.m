function [meshData,IEN,vertexData,IBC,BCval,totalNode,totalDOF] = ...
    createMesh(TotalRow,TotalColumn)
%TotalRow = 5;
%TotalColumn = 5;
global Grid_size;
domainx = [0,1];
domainy = [0,1];
nsd = 2;

BCuleft=0;BCuright=0;BCutop=1;BCubuttom=0;
BCvleft=0;BCvright=0;BCvtop=0;BCvbuttom=0;
BCpbuttom=0;

BCuV1=0;BCuV2=0;BCuV3=1;BCuV4=1;
BCvV1=0;BCvV2=0;BCvV3=0;BCvV4=0;
BCpV1=0;BCpV2=0;


Grid_size = (TotalRow-1)*(TotalColumn-1);
totalNode = (TotalRow)*(TotalColumn);
totalDOF = totalNode*(nsd+1);
meshData = zeros(Grid_size,4);
for i = 1:Grid_size
    [row,colume]=Index1To2(i,TotalRow-1);
    [index1] = Index2To1(row,colume,TotalRow);
    index2 = Index2To1(row,colume+1,TotalRow);
    index3 = Index2To1(row+1,colume+1,TotalRow);
    index4 = Index2To1(row+1,colume,TotalRow);
    meshData(i,:) = [index1,index2,index3,index4];
end
IEN = [meshData,meshData+totalNode,meshData+2*totalNode];

hx = (domainx(2)-domainx(1))/(TotalRow-1);
hy = (domainy(2)-domainy(1))/(TotalColumn-1);
vertexData = zeros(totalNode,2);
for i = 1:totalNode
    [row,colume]=Index1To2(i,TotalRow);
    vertexData(i,1) = (colume-1)*hx;
    vertexData(i,2) = (row-1)*hy;
end

%IBC
% top
IBC_u = [];
BCval_u = [];
IBC_v = [];
BCval_v=[];
IBC_p = [];
BCval_p=[];

[indexV1] = Index2To1(1,1,TotalRow);
[indexV2] = Index2To1(1,TotalColumn,TotalRow);
[indexV3] = Index2To1(TotalRow,TotalColumn,TotalRow);
[indexV4] = Index2To1(TotalRow,1,TotalRow);
IBC_u = [IBC_u,indexV1,indexV2,indexV3,indexV4];
BCval_u = [BCval_u,BCuV1,BCuV2,BCuV3,BCuV4];
IBC_v = [IBC_v,indexV1,indexV2,indexV3,indexV4];
BCval_v = [BCval_v,BCvV1,BCvV2,BCvV3,BCvV4];
IBC_p = [IBC_p,indexV1,indexV2];
BCval_p = [BCval_p,BCpV1,BCpV2];

rowN = TotalRow;
for colume = 2:TotalColumn-1
    [indexIBC] = Index2To1(rowN,colume,TotalRow);
    IBC_u = [IBC_u,indexIBC];
    BCval_u = [BCval_u,BCutop];
    IBC_v = [IBC_v,indexIBC];
    BCval_v = [BCval_v,BCvtop];
end
% buttom
rowN = 1;
for colume = 2:TotalColumn-1
    [indexIBC] = Index2To1(rowN,colume,TotalRow);
    IBC_u = [IBC_u,indexIBC];
    BCval_u = [BCval_u,BCubuttom];
    IBC_v = [IBC_v,indexIBC];
    BCval_v = [BCval_v,BCvbuttom];
    IBC_p = [IBC_p,indexIBC];
    BCval_p = [BCval_p,BCpbuttom];
end

%left
columeN = 1;
for row = 2:TotalRow-1
    [indexIBC] = Index2To1(row,columeN,TotalRow);
    IBC_u = [IBC_u,indexIBC];
    BCval_u = [BCval_u,BCuleft];
    IBC_v = [IBC_v,indexIBC];
    BCval_v = [BCval_v,BCvleft];
end
% right
columeN = TotalColumn;
for row = 2:TotalRow-1
    [indexIBC] = Index2To1(row,columeN,TotalRow);
    IBC_u = [IBC_u,indexIBC];
    BCval_u = [BCval_u,BCuright];
    IBC_v = [IBC_v,indexIBC];
    BCval_v = [BCval_v,BCvright];
end
IBC = [IBC_u,IBC_v+totalNode,IBC_p+2*totalNode];
BCval=[BCval_u,BCval_v,BCval_p];
end

function [row,column]=Index1To2(index2D,TotalRow)
row = floor((index2D-1)/TotalRow)+1;
column = index2D-(row-1)*TotalRow;
end

function [index2D] = Index2To1(row,colume,TotalRow)
index2D = (row-1)*TotalRow+colume;
end