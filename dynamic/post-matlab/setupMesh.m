function [pOrder,IENmesh,meshData,vertexData,IENall,pAll,...
    JwG,row,column,IBC_u,IBC_v,IBC_p,BC_u,BC_v,BC_p]=setupMesh(fileName)
global node_Size;
global TotalNode;
global Grid_size;
fileName
%% setup mesh
[Grid_size,node_Size,TotalNode,pOrder,maxDOF,mesh_ien,solution_ien,order_list,...
    vertexDataarray,BC_array,edgeNode,BCele] = readinMesh(fileName);
mesh = geo_Mesh(mesh_ien,vertexDataarray,BC_array,BCele,edgeNode,pOrder);
[meshData] = mesh.MeshData(); %mesh_ien
[vertexData] = mesh.vertexMesh(meshData); %vertex
[Area] = mesh.elementArea(meshData,vertexData); %element area
%[iper,IBC,BC_val,dBC_val,IBC_u,IBC_v,IBC_p,BC_u,BC_v,BC_p] = mesh.BoundaryCondition(vertexData); 
        
[BCB, lam_BCB,ndirec,ele_flux] = mesh.BoundaryFlux();
[BCF, lam_BCF, nFdirec] = mesh.BoundaryForce();

[IBC,IBC_u,IBC_v,IBC_p,BC_utag,BC_vtag,BC_ptag] = mesh.IBCHO();
[BC_val,BC_u,BC_v,BC_p] = mesh.BCvalHO(vertexData,IBC_u,IBC_v,IBC_p,...
            BC_utag,BC_vtag,BC_ptag,BCB,lam_BCB,mesh_ien,order_list);

    
n = nIntergerPoints(pOrder,0);
qPoints = TriGaussPoints(n);
nQ = size(qPoints,1);

IENmesh = zeros(Grid_size,3);
IENall = zeros(Grid_size,maxDOF);
pAll = zeros(Grid_size,3+1);
JwG = zeros(Grid_size,nQ);


for ele = 1:Grid_size 
    %if uniform p, use tabulated shapefunction and div values
    %else, calculate mixorder element shape functions     
    IEN_mesh = meshData{ele,2};
    ien = solution_ien(ele,:);
    IEN_all = ien(ien~=0);
    p_All = order_list(ele,:);

    IENmesh(ele,:) = IEN_mesh;
    IENall(ele,:) = IEN_all;
    pAll(ele,:) = p_All;
    
    qPoints = TriGaussPoints(n);    
    tri = SimplexGeometry(IEN_mesh,vertexData);
    [JInverse, detJ,gij,~,~,~] =  tri.jacobian();
    Jw = detJ*qPoints(:,3);
    JwG(ele,:) = Jw;
    
end
[row,column,Kindex,rowIndex,columnIndex,BCBindex] = K_index(IENall,IBC,BCB); 
end