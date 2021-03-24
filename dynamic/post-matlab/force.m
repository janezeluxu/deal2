function [Force] = force(order,num_start,nstp,mesh_file,bcele_file,ufile,pfile)
%% time stepping
tic
%% setup mesh
[solution_ien,p_ien,vertexData] = readin_mesh(mesh_file);
eleBC = readin_bcele(bcele_file);
%order = 1;
i = 0;
Force = zeros(nstp/10,2);
[BCF,lam_BCF] = getforceEle(eleBC);
for istp = num_start:nstp+num_start-1
    if mod(istp,10)==0
        i = i+1;        
        ufilename = strcat(ufile,'100-un-',string(istp+1),'.txt')
        pfilename = strcat(pfile,'100-pn-',string(istp+1),'.txt');
        uv = load(ufilename);
        p = load(pfilename);
        variable = [uv,p];
        [F] = calculateForce(variable,BCF,lam_BCF,solution_ien,p_ien,order,vertexData);
        Force(i,:) = F;
%         [Peclet,Pemin,Pemax] = getPecletNumber(u_n(1:TotalNode),u_n(TotalNode+1:2*TotalNode),Area,meshData,...
%             solution_ien,order_ien);
%         Pe(i,:) = [Pemin,Pemax];
%         %nonDtaum = gettaumNondimentional(u_n(1:TotalNode),u_n(TotalNode+1:2*TotalNode),tau_n,Area,meshData,...
%         %        solution_ien,order_ien);
%         statictau = getstaticTau(u_n(1:TotalNode),u_n(TotalNode+1:2*TotalNode),gijG,SF,sizea,sizeb,Area,meshData,...
%             solution_ien,order_ien);
    end
        
end
toc
end

function [BCF, lam] = getforceEle(eleBC)
BCF_edge = 50;
BoundaryELE = [];
lam = [];
for i = 1:size(eleBC,1)
    eID = eleBC(i,1);
    edge_tag = eleBC(i,2);
    
    for j = 1:size(BCF_edge,1)
        if (edge_tag == BCF_edge(j))
            BoundaryELE = [BoundaryELE,eID+1];
            lam = [lam, eleBC(i,3)];
        end
    end
end
BCF = BoundaryELE;
end