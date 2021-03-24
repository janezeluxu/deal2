function create_mesh(fileName,fileNameVertex,cellLevel,IENfilename,bcFile,outmeshfile,outbcfile)
clc
clear all
fileName = "../output/mesh/Shape1-1mesh.txt";
fileNameVertex = "../output/mesh/Shape1-1v_to_e_indices.txt";
cellLevel = load("prediction_Shape1_Re50");
IENfilename = "../output/mesh/Shape1-1EdgeInfo.txt";
bcFile = "../output/mesh/Shape1-1bc.txt";

outmeshfile = "../NS-steadyHO-quad/mesh/Shape1.txt";
outbcfile = "../NS-steadyHO-quad/mesh/Shape1bc.txt";
solutionFile = "../NS-steadyHO-quad/testvtk/cases/solution/cylinder/shape1.txt";
[solution_ien,p_ien,cellvertexDatax,cellvertexDatay,vertexData,v_e] ...
    = readin_mesh(fileName,fileNameVertex);
[Grid_size,vertices_size,Edge_size,total_DOF,meshIEN,mesh_Edge,Edge_v,IBC,IBC_edge] ...
    = load_mesh(IENfilename,bcFile);

uv = load("../output/solution/Shape1un-1-.txt");
p = load("../output/solution/Shape1pn-1-.txt");

v_count = vertices_size;
f_count = 1;
L_ev = ones(Edge_size,2)*2;
%cellLevel(1:10)
totalELE = 0;
vertexDatafinal = zeros(4000,2);
ufinal = zeros(4000,1);
vfinal = zeros(4000,1);
pfinal = zeros(4000,1);
for i = 1:Grid_size
    if cellLevel(i) >2
    cellLevel(i) = 2;
    end
end
%    cellLevel(10) = 1;
for i = 1:Grid_size
    n = cellLevel(i);
    totalELE = totalELE+2^(2*n);
end
totalELE
IEN = zeros(totalELE,4);
variable = [uv,p];

%for i = 1:Edge_size
for i = 1:Grid_size
    %i
    VertexList = meshIEN(i,:);
    %vID(1:4) = VertexList;
    EdgeList = mesh_Edge(i,:);
    nR = cellLevel(i);
    
    nEv = 2^nR+1;
    nV = (nEv)^2;
    vID = zeros(nV,1);
    v_coord = zeros(nV,1);
    for v = 1:length(VertexList)
        vID(v) = VertexList(v);
        v_coord(v,1) = vertexData(VertexList(v),1);
        v_coord(v,2) = vertexData(VertexList(v),2);
    end
    %vID
    eV = zeros(length(EdgeList),nEv);
    Lev = zeros(length(EdgeList),2);
    for e = 1:length(EdgeList)
        Lev(e,:) = L_ev(EdgeList(e),:);
        for j = 1:length(Edge_v(EdgeList(e),:))
            eV(e,j) = Edge_v(EdgeList(e),j);
        end
    end
    %nR
    
    %eV
    %% update vertex_info
    %nEv
%     if i == 1980
%         i
%         Lev
%         eV
%         nEv
%         vID
%         nR
%     end
    
    [vID] = updateExisting(Lev,eV,nEv,vID,nR);
%     if i == 154
%         vID
%         %vID(289),
%         %vid74 = vID(74),
%         %vID(4),vID(64)
%     end
    
    
    if nV>4
        for vv = 5:nV
            if vID(vv)==0
                vID(vv) = v_count+1;
                v_count = v_count+1;
            end
        end
    end
%     if i == 262
%         %vid289 = vID(289),
%         %vid74 = vID(74),
%         %vID(4),vID(64)
%     end
    %% update edge_info
    if nEv>2
        for edge = 1:4
            Lev1 = Lev(edge,1);
            
            if Lev1<nEv
                evstart = 5+(nEv-2)*(edge-1);
                eV(edge,3:nEv) = vID(evstart:evstart+nEv-3);
%                 if i == 1980
%                     edge
%                     evstart
%                     vID
%                     eV
%                 end
            end
    
        end
    end
    
%     if i == 1980
%         vID
%         eV
%     end

    for e = 1:length(EdgeList)
        if L_ev(EdgeList(e),1)>2
            L_ev(EdgeList(e),2) = nEv;
        else
            L_ev(EdgeList(e),1) = nEv;
        end
    end
    
    for e = 1:length(EdgeList)
        for j = 1:nEv
            Edge_v(EdgeList(e),j) = eV(e,j);
        end
    end
    
%     for e = 1:length(EdgeList)
%         if EdgeList(e) == 5277
%             i
%             eV
%         end
%     end
     %% update v_coord
    v_coordf = get_cordinates(v_coord,nEv,nR);
    %ele_v_coord(i,:) = v_coord;
    for vc = 1:nV
        vertexDatafinal(vID(vc),:) = v_coordf(vc,:);
    end
    %size(vertexData)
    
    %% update face_info
    IEN_ele = updateFace(vID,nR);
    nF = (nEv-1)^2;
    for face = 1:nF
        IEN(f_count,:) = IEN_ele(face,:);
        f_count = f_count+1;
    end

%     if i == 153
%         i
%         v_coordf
%         vID
%         Lev
%         eV
%         IEN_ele
%         f_count
%     end
%     if i == 2669
%         i
%         vID
%         Lev
%         eV
%         v_coordf
%         v3392 = vertexDatafinal(3392,:)
%         v3393 = vertexDatafinal(3393,:)
%     end
    
    
%     vID(10:12)
%     v_coord(10,:)
%     v_coord(11,:)
%     v_coord(12,:)
%     eV
%     IEN_ele
    %vertexData(3391,:)
    %vertexData(3392,:)
    
    [IENall,pAll] = elementIEN(i,solution_ien,p_ien);
    nssl = length(IENall);    
    variable_ele = zeros(nssl,1);
    for data = 1:nssl
        variable_ele(data) = variable(IENall(data));
    end
    [u,v,uv,p] = getuvp(variable_ele,pAll);
    uall = get_cordinates(u,nEv,nR);
    vall = get_cordinates(v,nEv,nR);
    pall = get_cordinates(p,nEv,nR);
    
    for vc = 1:nV
        ufinal(vID(vc),:) = uall(vc,:);
        vfinal(vID(vc),:) = vall(vc,:);
        pfinal(vID(vc),:) = pall(vc,:);
    end
    
end
%v3392 = vertexDatafinal(3392,:)
%v3393 = vertexDatafinal(3393,:)
%v3394 = vertexDatafinal(3394,:)
%IBC_edge
%L_ev
solution = [ufinal;vfinal;pfinal];
HanNode = [];
HanCoef = [];
for i = 1:Edge_size
    if L_ev(i,1)~= L_ev(i,2)
        %i
        %L_ev(i,:)
        %vID = Edge_v(i,:)
        [Node,Coef] = getHanNode(L_ev(i,1),L_ev(i,2),Edge_v(i,:));
        HanNode = [HanNode;Node];
        HanCoef = [HanCoef;Coef];
        
%         if Node(1) == 42369
%             i
%             L_ev(i,:)
%             vID = Edge_v(i,:)
%             Node
%             Coef
%         end
    end
end
v_count
%% write to file
%% IEN, vertexData, IBC,HanNode,HanCoef
fileID = fopen(outmeshfile,'w');
fprintf(fileID,'%6s %d\n','Grid_size',totalELE);
fprintf(fileID,'%6s %d\n','Node_size',v_count);
fprintf(fileID,'%6s %d\n','Total_DOF',v_count);
fprintf(fileID,'%6s %d\n','maxDOF',4);
fprintf(fileID,'%6s \n','Solution_Order');
fprintf(fileID,'%d \n',1);
fprintf(fileID,'%6s \n','meshIEN');
for f_count = 1:totalELE
    fprintf(fileID,'%d %d %d %d \n',IEN(f_count,:));
end
fprintf(fileID,'%6s \n','solutionIEN');
for f_count = 1:totalELE
    fprintf(fileID,'%d %d %d %d\n',IEN(f_count,:));
end
fprintf(fileID,'%6s \n','solutionOrder');
for f_count = 1:totalELE
    fprintf(fileID,'%d %d %d %d %d\n',1,1,1,1,1);
end
fprintf(fileID,'%6s \n','vertexCoord');
for vv = 1:v_count
    fprintf(fileID,'%2.10f %2.10f \n',vertexDatafinal(vv,:));
end
fclose(fileID);

fileID = fopen(outbcfile,'w');
fprintf(fileID,'%s %d \n','ibc',size(IBC,1));
%size(IBC,1)
for v_count = 1:size(IBC,1)
    fprintf(fileID,'%d %d \n',IBC(v_count,:));
end
%size(HanNode,1)
fprintf(fileID,'%s %d \n','HanNode',size(HanNode,1));
for v_count = 1:size(HanNode,1)
    fprintf(fileID,'%d %d %d\n',HanNode(v_count,:));
end
%size(HanCoef,1)
fprintf(fileID,'%s %d \n','HanCoef',size(HanCoef,1));
for v_count = 1:size(HanCoef,1)
    fprintf(fileID,'%2.10f %2.10f \n',HanCoef(v_count,:));
end
fclose(fileID);

fileID = fopen(solutionFile,'w');
formatSpec = '%2.16f \n';
fprintf(fileID,formatSpec,solution);
fclose(fileID);
end

function [NNode,Coef] = getHanNode(lEv1,lEv2,vID)
more = max(lEv1,lEv2);
less = min(lEv1,lEv2);

if less ==2 %% all nodes are relate to the end nodes
    nNode = more;
    Node = zeros(nNode-2,1);
    Coef1 = zeros(nNode-2,1);
    Node1 = zeros(nNode-2,1);
    Node2 = zeros(nNode-2,1);
    for i = 3:nNode
        Node(i-2) = vID(i);
        Coef1(i-2) = 1-(i-2)/(nNode-1);
        Node1(i-2) = vID(1);
        Node2(i-2) = vID(2);
    end
    Coef2 = 1-Coef1;
end

if less ==3 %% one size have 3 node
    if more ==5
        Node(1) = vID(3);
        Coef1(1) = 0.5;
        Coef2(1) = 0.5;
        Node1(1) = vID(1);
        Node2(1) = vID(4);
        
        Node(2) = vID(5);
        Coef1(2) = 0.5;
        Coef2(2) = 0.5;
        Node1(2) = vID(4);
        Node2(2) = vID(2);
        
        Node = Node';
        Node1 = Node1';
        Node2 = Node2';
        Coef1 = Coef1';
        Coef2 = Coef2';
    end
    
    if more ==9
        Node = zeros(6,1);
        Coef1 = zeros(6,1);
        Node1 = zeros(6,1);
        Node2 = zeros(6,1);
    
        nNode = 4;
        for i = 1:3
            Node(i) = vID(i+2);
            Node1(i) = vID(1);
            Node2(i) = vID(6);
            Coef1(i) = 1-(i)/nNode;
        end
        
        for i = 4:6
            Node(i) = vID(i+3);
            Node1(i) = vID(6);
            Node2(i) = vID(2);
            Coef1(i) = 1-(i-3)/nNode;
        end
        Coef2 = 1-Coef1;
    end
    
    if more ==17
        Node = zeros(14,1);
        Coef1 = zeros(14,1);
        Node1 = zeros(14,1);
        Node2 = zeros(14,1);
    
        nNode = 8;
        for i = 1:7
            Node(i) = vID(i+2);
            Node1(i) = vID(1);
            Node2(i) = vID(10);
            Coef1(i) = 1-(i)/nNode;
        end
        
        for i = 8:14
            Node(i) = vID(i+3);
            Node1(i) = vID(10);
            Node2(i) = vID(2);
            Coef1(i) = 1-(i-7)/nNode;
        end
        Coef2 = 1-Coef1;
        vID
        NNode = [Node,Node1,Node2];
        Coef = [Coef1,Coef2];
        
    end
end

if less ==5 %% one size have 5 node
    if more ==9
        for i = 1:4
            Node(i) = vID(2*i+1);
            Node1(i) = vID(2*i);
            Coef1(i) = 0.5;
        end
        Node1(1) = vID(1);
        Node2(4) = vID(2);
        for i = 1:3
            Node2(i) = vID(2*(i+1));
        end
        
        Coef2 = 1-Coef1;
        Node = Node';
        Node1 = Node1';
        Node2 = Node2';
        Coef1 = Coef1';
        Coef2 = Coef2';
    end
    
    if more ==17
        Node = [vID(3);vID(4);vID(5);vID(7);vID(8);vID(9);...
            vID(11);vID(12);vID(13);vID(15);vID(16);vID(17)];
        Node1 = [vID(1);vID(1);vID(1);vID(6);vID(6);vID(6);...
            vID(10);vID(10);vID(10);vID(14);vID(14);vID(14);];
        Node2 = [vID(6);vID(6);vID(6);vID(10);vID(10);vID(10);...
            vID(14);vID(14);vID(14);vID(2);vID(2);vID(2);];
        Coef1 = [0.75;0.5;0.25;0.75;0.5;0.25;0.75;0.5;0.25;0.75;0.5;0.25];
        Coef2 = 1-Coef1;
    end
    
end

if less ==9 %% one size have 9 node
    if more ==17
        for i = 1:8
            Node(i) = vID(2*i+1);
            Node1(i) = vID(2*i);
            Coef1(i) = 0.5;
        end
        Node1(1) = vID(1);
        Node2(8) = vID(2);
        for i = 1:7
            Node2(i) = vID(2*(i+1));
        end
        
        Coef2 = 1-Coef1;
        Node = Node';
        Node1 = Node1';
        Node2 = Node2';
        Coef1 = Coef1';
        Coef2 = Coef2';
    end  
end

NNode = [Node,Node1,Node2];
Coef = [Coef1,Coef2];
end

function [vID] = updateExisting(lEv_edge,eV,nEv,vID,level)
for e = 1:4
    if nEv>2
        lEv = lEv_edge(e,1);
        levlevel = log2(lEv-1);
        if nEv <=lEv %% all edge node are existing
            Reorder_eV = reorder_edge_node(eV(e,:),levlevel,nEv); %% order edge node to levels
            evstart = 5+(nEv-2)*(e-1);
            vID(evstart:evstart+nEv-3) = convert_back(Reorder_eV(1:nEv-2),level,nEv); %% convert back
        else %% some or all edge node will need to be created
            if lEv>2
                %levlevel = log2(lEv-1);
                Reorder_eV = reorder_edge_node(eV(e,:),levlevel,lEv); %% order edge node to levels
                evstart = 5+(nEv-2)*(e-1);
                temp = Reorder_eV(1:lEv-2);
                tempIndex = reorder_index(level)-3;
                for index = 1:lEv-2
                   vID(evstart+tempIndex(index)) = temp(index);
                end
            end
        end
    end
end
%vID
end

function Reorder_eV = reorder_edge_node(eV,level,lEv)
Reorder_index = reorder_index(level);
Reorder_eV = zeros(lEv-2,1);
for i = 1:lEv-2
    index = Reorder_index(i);
    Reorder_eV(i) = eV(index);
end
%Reorder_index
%Reorder_eV = Reorder_eV(3:end);
end

function converted_eV = convert_back(eV,level,lEv)
Reorder_index = reorder_index(level);
converted_eV = zeros(lEv-2,1);
for i = 1:lEv-2
    index = Reorder_index(i);
    converted_eV(index-2) = eV(i);
end
%converted_eV = converted_eV(3:end)
end

function  Reorder_index = reorder_index(level)
%% reorder edge node based on levels

temp = zeros(level,1);
for i = 2:level+1
    temp(i-1) = 2^(i-2);
end

s = 3;
start = zeros(level,1);
interval = zeros(level,1);
for i = 1:level
    s = s+temp(i);
    start(i+1) = s;
    interval(i) = 2^i; 
end
start(1) = 3;

if level>0
    count = 1;
    for i = level:-1:1
        j = 2^(i-1)+2;
        k = 2^i+1;
        
        ss = start(count);
        in = interval(count);
        
        c = 0;
        for ll = j:k
            Reorder_index(ll) = ss+in*c;
            c = c+1;
        end
        count = count+1;
    end
end

Reorder_index = Reorder_index(3:end);
end

function v_coord = get_cordinates(v_coord,nEv,nR)
%nEv
tmpindex = [1,3;2,4;1,2;3,4];
if nEv>2
    for edge = 1:4
        evstart = 5+(nEv-2)*(edge-1);
        %evstart:evstart+nEv-3
        cc = 1/(nEv-1);
        count = 1;
        index1 = tmpindex(edge,1);
        index2 = tmpindex(edge,2);
        for vv = evstart:evstart+nEv-3
            c2 = count*cc;
            c1 = 1-c2;
            v_coord(vv,:) = c1*v_coord(index1,:)+c2*v_coord(index2,:);
            %eV(edge,3:nEv) = vID(evstart:evstart+nEv-3);
            count = count+1;
        end
    end
    
    %nFace = (nEv-2)^2;
    if nR == 1
        v_coord(9,:) = 0.25*(v_coord(5,:)+v_coord(6,:)+v_coord(7,:)+v_coord(8,:));
    end
    
    if nR ==2
        v_coord(21,:) = 0.25*(v_coord(6,:)+v_coord(9,:)+v_coord(12,:)+v_coord(15,:));
        
        v_coord(18,:) = 0.5*(v_coord(12,:)+v_coord(21,:));
        v_coord(20,:) = 0.5*(v_coord(6,:)+v_coord(21,:));
        v_coord(22,:) = 0.5*(v_coord(21,:)+v_coord(9,:));
        v_coord(24,:) = 0.5*(v_coord(21,:)+v_coord(15,:));
        
        v_coord(17,:) = 0.25*(v_coord(5,:)+v_coord(18,:)+v_coord(11,:)+v_coord(20,:));
        v_coord(19,:) = 0.25*(v_coord(18,:)+v_coord(8,:)+v_coord(13,:)+v_coord(22,:));
        v_coord(23,:) = 0.25*(v_coord(20,:)+v_coord(14,:)+v_coord(7,:)+v_coord(24,:));
        v_coord(25,:) = 0.25*(v_coord(22,:)+v_coord(16,:)+v_coord(24,:)+v_coord(10,:));

    end
    
    if nR ==3
                v_coord(57,:) = 0.25*(v_coord(8,:)+v_coord(15,:)+v_coord(22,:)+v_coord(29,:));
        
        v_coord(55,:) = 0.5*(v_coord(8,:)+v_coord(57,:));
        v_coord(59,:) = 0.5*(v_coord(15,:)+v_coord(57,:));
        v_coord(43,:) = 0.5*(v_coord(22,:)+v_coord(57,:));
        v_coord(71,:) = 0.5*(v_coord(29,:)+v_coord(57,:));
        
        v_coord(41,:) = 0.25*(v_coord(6,:)+v_coord(43,:)+v_coord(20,:)+v_coord(55,:));
        v_coord(45,:) = 0.25*(v_coord(43,:)+v_coord(13,:)+v_coord(24,:)+v_coord(59,:));
        v_coord(69,:) = 0.25*(v_coord(10,:)+v_coord(71,:)+v_coord(55,:)+v_coord(27,:));
        v_coord(73,:) = 0.25*(v_coord(71,:)+v_coord(17,:)+v_coord(59,:)+v_coord(31,:));
        
         v_coord(34,:) = 0.5*(v_coord(20,:)+v_coord(41,:));
         v_coord(36,:) = 0.5*(v_coord(22,:)+v_coord(43,:));
         v_coord(38,:) = 0.5*(v_coord(24,:)+v_coord(45,:));  
         
         v_coord(40,:) = 0.5*(v_coord(6,:)+v_coord(41,:));
         v_coord(42,:) = 0.5*(v_coord(41,:)+v_coord(43,:));
         v_coord(44,:) = 0.5*(v_coord(43,:)+v_coord(45,:));
         v_coord(46,:) = 0.5*(v_coord(45,:)+v_coord(13,:));
         
         v_coord(48,:) = 0.5*(v_coord(41,:)+v_coord(55,:));
         v_coord(50,:) = 0.5*(v_coord(43,:)+v_coord(57,:));
         v_coord(52,:) = 0.5*(v_coord(45,:)+v_coord(59,:));  
        
         v_coord(54,:) = 0.5*(v_coord(8,:)+v_coord(55,:));
         v_coord(56,:) = 0.5*(v_coord(55,:)+v_coord(57,:));
         v_coord(58,:) = 0.5*(v_coord(57,:)+v_coord(59,:));
         v_coord(60,:) = 0.5*(v_coord(59,:)+v_coord(15,:));
         
         v_coord(62,:) = 0.5*(v_coord(55,:)+v_coord(69,:));
         v_coord(64,:) = 0.5*(v_coord(57,:)+v_coord(71,:));
         v_coord(66,:) = 0.5*(v_coord(59,:)+v_coord(73,:));  
        
         v_coord(68,:) = 0.5*(v_coord(10,:)+v_coord(69,:));       
         v_coord(70,:) = 0.5*(v_coord(69,:)+v_coord(71,:));
         v_coord(72,:) = 0.5*(v_coord(71,:)+v_coord(73,:));
         v_coord(74,:) = 0.5*(v_coord(73,:)+v_coord(17,:));
         
         v_coord(76,:) = 0.5*(v_coord(69,:)+v_coord(27,:));
         v_coord(78,:) = 0.5*(v_coord(71,:)+v_coord(29,:));
         v_coord(80,:) = 0.5*(v_coord(73,:)+v_coord(31,:));  
        
        v_coord(33,:) = 0.25*(v_coord(19,:)+v_coord(40,:)+v_coord(5,:)+v_coord(34,:));
        v_coord(35,:) = 0.25*(v_coord(21,:)+v_coord(42,:)+v_coord(34,:)+v_coord(36,:));
        v_coord(37,:) = 0.25*(v_coord(23,:)+v_coord(44,:)+v_coord(36,:)+v_coord(38,:));
        v_coord(39,:) = 0.25*(v_coord(25,:)+v_coord(46,:)+v_coord(38,:)+v_coord(12,:));

        v_coord(47,:) = 0.25*(v_coord(7,:)+v_coord(48,:)+v_coord(40,:)+v_coord(54,:));
        v_coord(49,:) = 0.25*(v_coord(48,:)+v_coord(50,:)+v_coord(42,:)+v_coord(56,:));
        v_coord(51,:) = 0.25*(v_coord(50,:)+v_coord(52,:)+v_coord(44,:)+v_coord(58,:));
        v_coord(53,:) = 0.25*(v_coord(52,:)+v_coord(14,:)+v_coord(46,:)+v_coord(60,:));
        
        v_coord(61,:) = 0.25*(v_coord(9,:)+v_coord(62,:)+v_coord(54,:)+v_coord(68,:));
        v_coord(63,:) = 0.25*(v_coord(70,:)+v_coord(56,:)+v_coord(62,:)+v_coord(64,:));
        v_coord(65,:) = 0.25*(v_coord(58,:)+v_coord(72,:)+v_coord(64,:)+v_coord(66,:));
        v_coord(67,:) = 0.25*(v_coord(60,:)+v_coord(74,:)+v_coord(66,:)+v_coord(16,:));
        
        v_coord(75,:) = 0.25*(v_coord(68,:)+v_coord(26,:)+v_coord(11,:)+v_coord(76,:));
        v_coord(77,:) = 0.25*(v_coord(28,:)+v_coord(70,:)+v_coord(76,:)+v_coord(78,:));
        v_coord(79,:) = 0.25*(v_coord(30,:)+v_coord(72,:)+v_coord(78,:)+v_coord(80,:));
        v_coord(81,:) = 0.25*(v_coord(32,:)+v_coord(74,:)+v_coord(80,:)+v_coord(18,:));
    end
    
    if nR ==4
        v_coord(177,:) = 0.25*(v_coord(12,:)+v_coord(27,:)+v_coord(42,:)+v_coord(57,:));
        
        v_coord(173,:) = 0.5*(v_coord(12,:)+v_coord(177,:));
        v_coord(181,:) = 0.5*(v_coord(27,:)+v_coord(177,:));
        v_coord(117,:) = 0.5*(v_coord(42,:)+v_coord(177,:));
        v_coord(237,:) = 0.5*(v_coord(57,:)+v_coord(177,:));
        
        v_coord(113,:) = 0.25*(v_coord(8,:)+v_coord(117,:)+v_coord(38,:)+v_coord(173,:));
        v_coord(121,:) = 0.25*(v_coord(117,:)+v_coord(23,:)+v_coord(46,:)+v_coord(181,:));
        v_coord(233,:) = 0.25*(v_coord(16,:)+v_coord(237,:)+v_coord(173,:)+v_coord(53,:));
        v_coord(241,:) = 0.25*(v_coord(237,:)+v_coord(31,:)+v_coord(181,:)+v_coord(61,:));
        
        v_coord(83,:) = 0.5*(v_coord(38,:)+v_coord(113,:));
        v_coord(87,:) = 0.5*(v_coord(42,:)+v_coord(117,:));
        v_coord(91,:) = 0.5*(v_coord(46,:)+v_coord(121,:));  
        
        v_coord(111,:) = 0.5*(v_coord(8,:)+v_coord(113,:));
        v_coord(115,:) = 0.5*(v_coord(113,:)+v_coord(117,:));
        v_coord(119,:) = 0.5*(v_coord(117,:)+v_coord(121,:));
        v_coord(123,:) = 0.5*(v_coord(121,:)+v_coord(23,:));
        
        v_coord(143,:) = 0.5*(v_coord(113,:)+v_coord(173,:));
        v_coord(147,:) = 0.5*(v_coord(117,:)+v_coord(177,:));
        v_coord(151,:) = 0.5*(v_coord(121,:)+v_coord(181,:));  
        
        v_coord(171,:) = 0.5*(v_coord(12,:)+v_coord(173,:));
        v_coord(175,:) = 0.5*(v_coord(173,:)+v_coord(177,:));
        v_coord(179,:) = 0.5*(v_coord(177,:)+v_coord(181,:));
        v_coord(183,:) = 0.5*(v_coord(181,:)+v_coord(27,:));
        
        v_coord(203,:) = 0.5*(v_coord(173,:)+v_coord(233,:));
        v_coord(207,:) = 0.5*(v_coord(177,:)+v_coord(237,:));
        v_coord(211,:) = 0.5*(v_coord(181,:)+v_coord(241,:));  
        
        v_coord(231,:) = 0.5*(v_coord(16,:)+v_coord(233,:));       
        v_coord(235,:) = 0.5*(v_coord(233,:)+v_coord(237,:));
        v_coord(239,:) = 0.5*(v_coord(237,:)+v_coord(241,:));
        v_coord(243,:) = 0.5*(v_coord(241,:)+v_coord(31,:));
        
        v_coord(263,:) = 0.5*(v_coord(233,:)+v_coord(53,:));
        v_coord(267,:) = 0.5*(v_coord(237,:)+v_coord(57,:));
        v_coord(271,:) = 0.5*(v_coord(241,:)+v_coord(61,:));  
        
        v_coord(81,:) = 0.25*(v_coord(6,:)+v_coord(83,:)+v_coord(36,:)+v_coord(111,:));
        v_coord(85,:) = 0.25*(v_coord(83,:)+v_coord(87,:)+v_coord(40,:)+v_coord(115,:));
        v_coord(89,:) = 0.25*(v_coord(87,:)+v_coord(91,:)+v_coord(119,:)+v_coord(44,:));
        v_coord(93,:) = 0.25*(v_coord(91,:)+v_coord(21,:)+v_coord(48,:)+v_coord(123,:));

        v_coord(141,:) = 0.25*(v_coord(10,:)+v_coord(143,:)+v_coord(111,:)+v_coord(171,:));
        v_coord(145,:) = 0.25*(v_coord(143,:)+v_coord(147,:)+v_coord(115,:)+v_coord(175,:));
        v_coord(149,:) = 0.25*(v_coord(147,:)+v_coord(151,:)+v_coord(119,:)+v_coord(179,:));
        v_coord(153,:) = 0.25*(v_coord(151,:)+v_coord(25,:)+v_coord(123,:)+v_coord(183,:));
        
        v_coord(201,:) = 0.25*(v_coord(14,:)+v_coord(203,:)+v_coord(171,:)+v_coord(231,:));
        v_coord(205,:) = 0.25*(v_coord(203,:)+v_coord(207,:)+v_coord(175,:)+v_coord(235,:));
        v_coord(209,:) = 0.25*(v_coord(207,:)+v_coord(211,:)+v_coord(179,:)+v_coord(239,:));
        v_coord(213,:) = 0.25*(v_coord(211,:)+v_coord(29,:)+v_coord(183,:)+v_coord(243,:));
        
        v_coord(261,:) = 0.25*(v_coord(18,:)+v_coord(263,:)+v_coord(231,:)+v_coord(51,:));
        v_coord(265,:) = 0.25*(v_coord(263,:)+v_coord(267,:)+v_coord(235,:)+v_coord(55,:));
        v_coord(269,:) = 0.25*(v_coord(267,:)+v_coord(271,:)+v_coord(239,:)+v_coord(59,:));
        v_coord(273,:) = 0.25*(v_coord(271,:)+v_coord(33,:)+v_coord(243,:)+v_coord(63,:));
        
        for i = 0:6
            v_coord(66+2*i,:) = 0.5*v_coord(36+2*i,:)+0.5*v_coord(81+2*i,:);
            v_coord(96+2*i,:) = 0.5*v_coord(81+2*i,:)+0.5*v_coord(111+2*i,:);
            v_coord(126+2*i,:) = 0.5*v_coord(111+2*i,:)+0.5*v_coord(141+2*i,:);
            v_coord(156+2*i,:) = 0.5*v_coord(141+2*i,:)+0.5*v_coord(171+2*i,:);
            v_coord(186+2*i,:) = 0.5*v_coord(171+2*i,:)+0.5*v_coord(201+2*i,:);
            v_coord(216+2*i,:) = 0.5*v_coord(201+2*i,:)+0.5*v_coord(231+2*i,:);
            v_coord(246+2*i,:) = 0.5*v_coord(231+2*i,:)+0.5*v_coord(261+2*i,:);
            v_coord(276+2*i,:) = 0.5*v_coord(261+2*i,:)+0.5*v_coord(51+2*i,:);
        end
        
        v_coord(80,:) = 0.5*v_coord(6,:)+0.5*v_coord(81,:);
        v_coord(110,:) = 0.5*v_coord(8,:)+0.5*v_coord(111,:);
        v_coord(140,:) = 0.5*v_coord(10,:)+0.5*v_coord(141,:);
        v_coord(170,:) = 0.5*v_coord(12,:)+0.5*v_coord(171,:);
        v_coord(200,:) = 0.5*v_coord(14,:)+0.5*v_coord(201,:);
        v_coord(230,:) = 0.5*v_coord(16,:)+0.5*v_coord(231,:);
        v_coord(260,:) = 0.5*v_coord(18,:)+0.5*v_coord(261,:);
        for i = 0:5
            v_coord(82+2*i,:) = 0.5*v_coord(81+2*i,:)+0.5*v_coord(83+2*i,:);
            v_coord(112+2*i,:) = 0.5*v_coord(111+2*i,:)+0.5*v_coord(113+2*i,:);
            v_coord(142+2*i,:) = 0.5*v_coord(141+2*i,:)+0.5*v_coord(143+2*i,:);
            v_coord(172+2*i,:) = 0.5*v_coord(171+2*i,:)+0.5*v_coord(173+2*i,:);
            v_coord(202+2*i,:) = 0.5*v_coord(201+2*i,:)+0.5*v_coord(203+2*i,:);
            v_coord(232+2*i,:) = 0.5*v_coord(231+2*i,:)+0.5*v_coord(233+2*i,:);
            v_coord(262+2*i,:) = 0.5*v_coord(261+2*i,:)+0.5*v_coord(263+2*i,:);
        end
        v_coord(94,:) = 0.5*v_coord(93,:)+0.5*v_coord(21,:);
        v_coord(124,:) = 0.5*v_coord(123,:)+0.5*v_coord(23,:);
        v_coord(154,:) = 0.5*v_coord(153,:)+0.5*v_coord(25,:);
        v_coord(184,:) = 0.5*v_coord(183,:)+0.5*v_coord(27,:);
        v_coord(214,:) = 0.5*v_coord(213,:)+0.5*v_coord(29,:);
        v_coord(244,:) = 0.5*v_coord(243,:)+0.5*v_coord(31,:);
        v_coord(274,:) = 0.5*v_coord(273,:)+0.5*v_coord(33,:);

        v_coord(65,:) = 0.25*(v_coord(5,:)+v_coord(66,:)+v_coord(35,:)+v_coord(80,:));
        v_coord(95,:) = 0.25*(v_coord(7,:)+v_coord(96,:)+v_coord(80,:)+v_coord(110,:));
        v_coord(125,:) = 0.25*(v_coord(9,:)+v_coord(126,:)+v_coord(110,:)+v_coord(140,:));
        v_coord(155,:) = 0.25*(v_coord(11,:)+v_coord(156,:)+v_coord(140,:)+v_coord(170,:));
        v_coord(185,:) = 0.25*(v_coord(13,:)+v_coord(186,:)+v_coord(170,:)+v_coord(200,:));
        v_coord(215,:) = 0.25*(v_coord(15,:)+v_coord(216,:)+v_coord(200,:)+v_coord(230,:));
        v_coord(245,:) = 0.25*(v_coord(17,:)+v_coord(246,:)+v_coord(230,:)+v_coord(260,:));
        v_coord(275,:) = 0.25*(v_coord(19,:)+v_coord(276,:)+v_coord(260,:)+v_coord(50,:));
        for i = 0:5
            v_coord(67+2*i,:) = 0.25*(v_coord(66+2*i,:)+v_coord(68+2*i,:)+v_coord(37+2*i,:)+v_coord(82+2*i,:));
            v_coord(97+2*i,:) = 0.25*(v_coord(96+2*i,:)+v_coord(98+2*i,:)+v_coord(82+2*i,:)+v_coord(112+2*i,:));
            v_coord(127+2*i,:) = 0.25*(v_coord(126+2*i,:)+v_coord(128+2*i,:)+v_coord(112+2*i,:)+v_coord(142+2*i,:));
            v_coord(157+2*i,:) = 0.25*(v_coord(156+2*i,:)+v_coord(158+2*i,:)+v_coord(142+2*i,:)+v_coord(172+2*i,:));
            v_coord(187+2*i,:) = 0.25*(v_coord(186+2*i,:)+v_coord(188+2*i,:)+v_coord(172+2*i,:)+v_coord(202+2*i,:));
            v_coord(217+2*i,:) = 0.25*(v_coord(216+2*i,:)+v_coord(218+2*i,:)+v_coord(202+2*i,:)+v_coord(232+2*i,:));
            v_coord(247+2*i,:) = 0.25*(v_coord(246+2*i,:)+v_coord(248+2*i,:)+v_coord(232+2*i,:)+v_coord(262+2*i,:));
            v_coord(277+2*i,:) = 0.25*(v_coord(276+2*i,:)+v_coord(278+2*i,:)+v_coord(262+2*i,:)+v_coord(52+2*i,:));
        end
        v_coord(79,:) = 0.25*(v_coord(78,:)+v_coord(20,:)+v_coord(49,:)+v_coord(94,:));
        v_coord(109,:) = 0.25*(v_coord(108,:)+v_coord(22,:)+v_coord(94,:)+v_coord(124,:));
        v_coord(139,:) = 0.25*(v_coord(138,:)+v_coord(24,:)+v_coord(124,:)+v_coord(154,:));
        v_coord(169,:) = 0.25*(v_coord(168,:)+v_coord(26,:)+v_coord(154,:)+v_coord(184,:));
        v_coord(199,:) = 0.25*(v_coord(198,:)+v_coord(28,:)+v_coord(184,:)+v_coord(214,:));
        v_coord(229,:) = 0.25*(v_coord(228,:)+v_coord(30,:)+v_coord(214,:)+v_coord(244,:));
        v_coord(259,:) = 0.25*(v_coord(258,:)+v_coord(32,:)+v_coord(244,:)+v_coord(274,:));
        v_coord(289,:) = 0.25*(v_coord(288,:)+v_coord(34,:)+v_coord(274,:)+v_coord(64,:));
    end
    
end
  
%v_coord
%plot(v_coord(:,1),v_coord(:,2),'bo')
end

function IEN = updateFace(vID,nR)
if nR == 0
    IEN = [vID(1),vID(2),vID(4),vID(3)];
end

if nR == 1
    IEN(1,:) = [vID(1),vID(7),vID(9),vID(5)];
    IEN(2,:) = [vID(7),vID(2),vID(6),vID(9)];
    IEN(3,:) = [vID(5),vID(9),vID(8),vID(3)];
    IEN(4,:) = [vID(9),vID(6),vID(4),vID(8)];
end

if nR == 2
    IEN(1,:) = [vID(1),vID(11),vID(17),vID(5)];
    IEN(2,:) = [vID(11),vID(12),vID(18),vID(17)];
    IEN(3,:) = [vID(12),vID(13),vID(19),vID(18)];
    IEN(4,:) = [vID(13),vID(2),vID(8),vID(19)];
    
    IEN(5,:) = [vID(5),vID(17),vID(20),vID(6)];
    IEN(6,:) = [vID(17),vID(18),vID(21),vID(20)];
    IEN(7,:) = [vID(18),vID(19),vID(22),vID(21)];
    IEN(8,:) = [vID(19),vID(8),vID(9),vID(22)];
    
    IEN(9,:) = [vID(6),vID(20),vID(23),vID(7)];
    IEN(10,:) = [vID(20),vID(21),vID(24),vID(23)];
    IEN(11,:) = [vID(21),vID(22),vID(25),vID(24)];
    IEN(12,:) = [vID(22),vID(9),vID(10),vID(25)];
    
    IEN(13,:) = [vID(7),vID(23),vID(14),vID(3)];
    IEN(14,:) = [vID(23),vID(24),vID(15),vID(14)];
    IEN(15,:) = [vID(24),vID(25),vID(16),vID(15)];
    IEN(16,:) = [vID(25),vID(10),vID(4),vID(16)];
end

if nR == 3
    IEN(1,:) = [vID(1),vID(19),vID(33),vID(5)];
    IEN(8,:) = [vID(25),vID(2),vID(12),vID(39)];
    for i = 0:5
        IEN(9+8*i,:) = [vID(5+i),vID(33+7*i),vID(40+7*i),vID(6+i)];
        IEN(16+8*i,:) = [vID(39+7*i),vID(12+i),vID(13+i),vID(46+7*i)];
        
        for j = 0:5
            IEN(10+i+8*j,:) = [vID(33+i+7*j),vID(34+i+7*j),vID(41+i+7*j),vID(40+i+7*j)];
        end
        
        IEN(2+i,:) = [vID(19+i),vID(20+i),vID(34+i),vID(33+i)];
        
%         IEN(10+i,:) = [vID(33+i),vID(34+i),vID(41+i),vID(40+i)];
%         IEN(18+i,:) = [vID(40+i),vID(41+i),vID(48+i),vID(47+i)];
%         IEN(26+i,:) = [vID(47+i),vID(48+i),vID(55+i),vID(54+i)];
%         IEN(34+i,:) = [vID(54+i),vID(55+i),vID(62+i),vID(61+i)];
%         IEN(42+i,:) = [vID(61+i),vID(62+i),vID(69+i),vID(68+i)];
%         IEN(50+i,:) = [vID(68+i),vID(69+i),vID(76+i),vID(75+i)];
        
        IEN(58+i,:) = [vID(75+i),vID(76+i),vID(27+i),vID(26+i)];        
    end
    IEN(57,:) = [vID(11),vID(75),vID(26),vID(3)];
    IEN(64,:) = [vID(81),vID(18),vID(4),vID(32)];

end

if nR == 4

    IEN(1,:) = [vID(1),vID(35),vID(65),vID(5)];
    IEN(16,:) = [vID(49),vID(2),vID(20),vID(79)];
    for i = 0:13
        IEN(17+16*i,:) = [vID(5+i),vID(65+15*i),vID(80+15*i),vID(6+i)];
        IEN(32+16*i,:) = [vID(79+15*i),vID(20+i),vID(21+i),vID(94+15*i)];
        IEN(2+i,:) = [vID(35+i),vID(36+i),vID(66+i),vID(65+i)];
        
        for j = 0:13
            IEN(18+i+16*j,:) = [vID(65+i+15*j),vID(66+i+15*j),vID(81+i+15*j),vID(80+i+15*j)];
        end
        
        IEN(242+i,:) = [vID(275+i),vID(276+i),vID(51+i),vID(50+i)];
    end
    IEN(241,:) = [vID(19),vID(275),vID(50),vID(3)];
    IEN(256,:) = [vID(289),vID(34),vID(4),vID(64)];
    
end

end

function [u,v,uv,p] = getuvp(variable_ele,pAll)
uvsize = 3*(pAll+1)^2;
if pAll <3
    u = variable_ele(1:3:uvsize);
    v = variable_ele(2:3:uvsize);
    p = variable_ele(3:3:end);
    uv = variable_ele([1:3:uvsize,2:3:uvsize]);
elseif pAll==3
    %1:3:uvsize
    u = variable_ele([1,4,7,10,13,14,19,20,25,26,31,32,37,38,39,40]);
    v = variable_ele([2,5,8,11,15,16,21,22,27,28,33,34,41,42,43,44]);
    p = variable_ele([3,6,9,12,17,18,23,24,29,30,35,36,45,46,47,48]);
    uv = variable_ele([1,4,7,10,13,14,19,20,25,26,31,32,37,38,39,40,...
        2,5,8,11,15,16,21,22,27,28,33,34,41,42,43,44]);
elseif pAll==5
    %1:3:uvsize
    u = variable_ele([1,4,7,10,13:16,25:28,37:40,49:52,61:76]);
    v = variable_ele([2,5,8,11,17:20,29:32,41:44,53:56,77:92]);
    p = variable_ele([3,6,9,12,21:24,33:36,45:48,57:60,93:108]);
    uv = variable_ele([1,4,7,10,13:16,25:28,37:40,49:52,61:76,...
        2,5,8,11,17:20,29:32,41:44,53:56,77:92]);
 elseif pAll==7
    u = variable_ele([1,4,7,10,13:18,31:36,49:54,67:72,85:120]);
    v = variable_ele([2,5,8,11,19:24,37:42,55:60,73:78,121:156]);
    p = variable_ele([3,6,9,12,25:30,43:48,61:66,79:84,157:192]);
    uv = variable_ele([1,4,7,10,13:18,31:36,49:54,67:72,85:120,...
        2,5,8,11,19:24,37:42,55:60,73:78,121:156]);   
end
end

function [IENall,pAll] = elementIEN(ele,solution_ien,p_ien)
%get the IEN array and P list for a given element
%IEN_mesh = ien_map(ele,:);
IENall = solution_ien(ele,:);
IENall = IENall(IENall~=0);
pAll = p_ien(ele);
%x = vertexData{1}(ele,:);
%y = vertexData{2}(ele,:);
%xy(1,:) = x;
%xy(2,:) = y;
end