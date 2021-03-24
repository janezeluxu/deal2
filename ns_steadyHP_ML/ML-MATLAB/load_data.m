function [input_feature,Grid_size] = ...
    load_data(fileName,fileNameVertex,uv,p)
%vertex_size =625;
%Grid_size = 576;
%p1 = 1;
%fileName = "../output/mesh/lidcavity-1mesh.txt";
%fileNameVertex = "../output/mesh/lidcavity-1v_to_e_indices.txt";
[IEN_all,p_All,cellvertexDatax,cellvertexDatay,vertexData,v_e] ...
    = readin_mesh(fileName,fileNameVertex);

Grid_size = length(p_All);


%uv = load("../output/solution/lidcavityRe400un-1-.txt");
%p = load("../output/solution/lidcavityRe400pn-1-.txt");
%cellType = load("cellType.txt");
variable = [uv,p];
input_feature = zeros(Grid_size,12);
for ele = 1:Grid_size
    %ele
    %% get p1 grid information, shape function evaluate at p2 grid quadrature points
    [IENall,pAll] = elementIEN(ele,IEN_all,p_All);
    nssl = length(IENall);    
    variable_ele = zeros(nssl,1);
    for i = 1:nssl
        variable_ele(i) = variable(IENall(i));
    end
    [u,v,uv,p] = getuvp(variable_ele,pAll);
    input_feature(ele,:) = [u',v',p'];
    %elecellType = cellType(ele)
end
% input_feature
% cellType
% mdl = fitcensemble(input_feature,cellType);
% mdl.resubLoss
% confusionchart(cellType,mdl.resubPredict)
% prediction = mdl.resubPredict;
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