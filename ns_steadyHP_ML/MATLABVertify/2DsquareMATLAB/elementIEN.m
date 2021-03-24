function[eleIEN,xcord,ycord,vIDs] = elementIEN(ele,meshData,vertexData,IEN)
vIDs = meshData(ele,:);
x1 = vertexData(vIDs(1),1);
y1 = vertexData(vIDs(1),2);
x2 = vertexData(vIDs(2),1);
y2 = vertexData(vIDs(2),2);
x3 = vertexData(vIDs(3),1);
y3 = vertexData(vIDs(3),2);
x4 = vertexData(vIDs(4),1);
y4 = vertexData(vIDs(4),2);
xcord = [x1,x2,x3,x4];
ycord = [y1,y2,y3,y4];
   
eleIEN = IEN(ele,:);
end

