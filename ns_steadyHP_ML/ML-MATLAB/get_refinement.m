function get_refinement()
clc
clear all

%cellLevel = load("prediction_Shape1_Re50");

cellLevel = load("ML_prediction_Shape1_Re50");
%ML_alg_lidcavity_Re800
EdgeFilename1 = "../output/mesh/Shape1-1EdgeInfo.txt";
bcFile = "../output/mesh/Shape1-1bc.txt";

vtkfile = '../output/paraview/Shape1-01.vtk';
fid = fopen(vtkfile, 'at');
ele_sizeo = length(cellLevel);
if fid ~= -1
    %fprintf(fid, ...... whatever.....
    fprintf ( fid, '\n CELL_DATA %d\n', ele_sizeo );  
    fprintf ( fid, '\n SCALARS levelarray int\n' );
    fprintf ( fid, 'LOOKUP_TABLE default\n' );
    
    for element = 1 : ele_sizeo
        fprintf ( fid, '  %d', cellLevel(element) );
    end
    fclose(fid);
end

[mesh_Edge,edge_face,IBC_edge] = getedgeinfo(EdgeFilename1,bcFile);
level_array = cellLevel;
refine_level = 1;
[ele_list,refine_flag]=get_refine_ele(level_array, refine_level,mesh_Edge,edge_face,IBC_edge);
filename = "../input/refine_element/Shape1-1cycle1.txt";
fileID = fopen(filename,'w');
for ele_count = 1:length(refine_flag)
    fprintf(fileID,'%d\n',refine_flag(ele_count));
end

last_cycle_level_array = level_array;
Refineflagname = "../output/mesh/refineFlag1.txt";
RefineFlag_temp = load(Refineflagname);
%size(RefineFlag_temp)
RefineFlag = RefineFlag_temp(1:2:end);
%size(RefineFlag)
non_refine_count=nnz(~RefineFlag);

%last_cycle_level_arraysize(last_cycle_level_array)
[level_array] = update_ele_level(last_cycle_level_array,...
   RefineFlag,non_refine_count);
sizelevelarray = size(level_array)

vtkfile = '../output/paraview/Shape1-02.vtk';
%vtkfile = '../output/paraview/lidcavityRe800-01.vtk';
ele_sizeo = length(level_array);
fid = fopen(vtkfile, 'at');
if fid ~= -1
    %fprintf(fid, ...... whatever.....
    fprintf ( fid, '\n CELL_DATA %d\n', ele_sizeo );  
    fprintf ( fid, '\n SCALARS level1 int\n' );
    fprintf ( fid, 'LOOKUP_TABLE default\n' );
    
    for element = 1 : ele_sizeo
        fprintf ( fid, '  %d', level_array(element) );
    end
    fclose(fid);
end

%% cycle 2
 EdgeFilename2 = "../output/mesh/Shape1-2EdgeInfo.txt";
 bcFile = "../output/mesh/Shape1-2bc.txt";
 [mesh_Edge,edge_face,IBC_edge] = getedgeinfo(EdgeFilename2,bcFile);
  refine_level = 2;
  [ele_list,refine_flag]=get_refine_ele(level_array, refine_level,mesh_Edge,edge_face,IBC_edge);
  filename = "../input/refine_element/Shape1-1cycle2.txt";
fileID = fopen(filename,'w');
for ele_count = 1:length(refine_flag)
    fprintf(fileID,'%d\n',refine_flag(ele_count));
end

Refineflagname = "../output/mesh/refineFlag2.txt";
RefineFlag_temp = load(Refineflagname);
RefineFlag = RefineFlag_temp(1:2:end);
non_refine_count=nnz(~RefineFlag);

last_cycle_level_array = level_array;
 [level_array] = update_ele_level(last_cycle_level_array,...
     RefineFlag,non_refine_count);
 
 vtkfile = '../output/paraview/Shape1-03.vtk';
ele_sizeo = length(level_array);
fid = fopen(vtkfile, 'at');
if fid ~= -1
    %fprintf(fid, ...... whatever.....
    fprintf ( fid, '\n CELL_DATA %d\n', ele_sizeo );  
    fprintf ( fid, '\n SCALARS level1 int\n' );
    fprintf ( fid, 'LOOKUP_TABLE default\n' );
    
    for element = 1 : ele_sizeo
        fprintf ( fid, '  %d', level_array(element) );
    end
    fclose(fid);
end

%% cycle 3
 EdgeFilename = "../output/mesh/Shape1-3EdgeInfo.txt";
 bcFile = "../output/mesh/Shape1-3bc.txt";
 [mesh_Edge,edge_face,IBC_edge] = getedgeinfo(EdgeFilename,bcFile);
  refine_level = 3;
  [ele_list,refine_flag]=get_refine_ele(level_array, refine_level,mesh_Edge,edge_face,IBC_edge);
  filename = "../input/refine_element/Shape1-1cycle3.txt";
fileID = fopen(filename,'w');
for ele_count = 1:length(refine_flag)
    fprintf(fileID,'%d\n',refine_flag(ele_count));
end

Refineflagname = "../output/mesh/refineFlag3.txt";
RefineFlag_temp = load(Refineflagname);
RefineFlag = RefineFlag_temp(1:2:end);
non_refine_count=nnz(~RefineFlag);

last_cycle_level_array = level_array;
 [level_array] = update_ele_level(last_cycle_level_array,...
     RefineFlag,non_refine_count);
 
 vtkfile = '../output/paraview/Shape1-04.vtk';
ele_sizeo = length(level_array);
fid = fopen(vtkfile, 'at');
if fid ~= -1
    %fprintf(fid, ...... whatever.....
    fprintf ( fid, '\n CELL_DATA %d\n', ele_sizeo );  
    fprintf ( fid, '\n SCALARS level1 int\n' );
    fprintf ( fid, 'LOOKUP_TABLE default\n' );
    
    for element = 1 : ele_sizeo
        fprintf ( fid, '  %d', level_array(element) );
    end
    fclose(fid);
end

%% cycle 4

 EdgeFilename = "../output/mesh/Shape1-4EdgeInfo.txt";
 bcFile = "../output/mesh/Shape1-4bc.txt";
 [mesh_Edge,edge_face,IBC_edge] = getedgeinfo(EdgeFilename,bcFile);
  refine_level = 4;
  [ele_list,refine_flag]=get_refine_ele(level_array, refine_level,mesh_Edge,edge_face,IBC_edge);
  filename = "../input/refine_element/Shape1-1cycle4.txt";
fileID = fopen(filename,'w');
for ele_count = 1:length(refine_flag)
    fprintf(fileID,'%d\n',refine_flag(ele_count));
end

Refineflagname = "../output/mesh/refineFlag4.txt";
RefineFlag_temp = load(Refineflagname);
RefineFlag = RefineFlag_temp(1:2:end);
non_refine_count=nnz(~RefineFlag);

last_cycle_level_array = level_array;
 [level_array] = update_ele_level(last_cycle_level_array,...
     RefineFlag,non_refine_count);
 
vtkfile = '../output/paraview/Shape1-05.vtk';
ele_sizeo = length(level_array);
fid = fopen(vtkfile, 'at');
if fid ~= -1
    %fprintf(fid, ...... whatever.....
    fprintf ( fid, '\n CELL_DATA %d\n', ele_sizeo );  
    fprintf ( fid, '\n SCALARS level1 int\n' );
    fprintf ( fid, 'LOOKUP_TABLE default\n' );
    
    for element = 1 : ele_sizeo
        fprintf ( fid, '  %d', level_array(element) );
    end
    fclose(fid);
end

end

function [ele_list,refine_flag]=get_refine_ele(level_array, ...
    refine_level,mesh_Edge,edge_face,IBC_edge)
nEle = length(level_array);
non_refine_count = 0;
ele_list = [];
refine_flag = zeros(nEle,1);
for i = 1:nEle
    edge_list = mesh_Edge(i,:);
    face_list = zeros(4,2);
    for e = 1:length(edge_list)
        edge = edge_list(e);
        face_list(e,:) = edge_face(edge,:);
    end
    [if_refine] = test_if_refine(i,level_array(i),refine_level,edge_list,face_list,IBC_edge);
    if (if_refine)
        ele_list = [ele_list,i];
        refine_flag(i) = 1;
    else
        non_refine_count = non_refine_count+1;
        refine_flag(i) = 0;
    end
end
refien_size = length(ele_list)
non_refine_count
refien_size*4+non_refine_count
end

function [level_array] = update_ele_level(last_cycle_level_array,RefineFlag,non_refine_count)
nEle = length(last_cycle_level_array)
ele_array = zeros(nEle,4);
%ele_edge_array = zeros(nEle,4);
%snon_refine_count = 1312
refien_size = nEle-non_refine_count;
total_next_level_ele = refien_size*4+non_refine_count;
count = 1;
%new_edge_count = max(max(last_cycle_mesh_Edge));
size(RefineFlag)
for ele  = 1:nEle
    %level = last_cycle_level_array(ele);
    if RefineFlag(ele)>0
        ele_array(ele,:) = [non_refine_count+1:non_refine_count+4];
        non_refine_count = non_refine_count+4;
        %ele_edge_array(ele,:) = [new_edge_count+1:new_edge_count+4];
        %new_edge_count = new_edge_count+4;
    else
        ele_array(ele,:) = [count,count,count,count];
        count = count+1;
    end
        
end
%ele_array'
non_refine_count
count
%ele_array'

level_array = zeros(total_next_level_ele,1);
%mesh_Edge = zeros(total_next_level_ele,4);
for i = 1:nEle
    level = last_cycle_level_array(i);
    %last_edge = last_cycle_mesh_Edge(i,:)
    %ele_edge = ele_edge_array(i,:)
    for j = 1:4
        ele = ele_array(i,j); 
        level_array(ele) = level;
       %mesh_Edge(ele,:) = get_ele_edge(last_edge,ele_edge);
    end   
end
%level_arr_size = size(level_array)
end



function [if_refine] = test_if_refine(ele,current_level, refine_level,edge_list,face_list,IBC_edge)
if current_level<refine_level
    if_refine = false;
else %% dont refine if it will create more than 1 non-conforming node
%     for i = 1:4
%         for j = 1:2
%             if face_list(i,j) == 0
%                 find(~face_list)
%                 edge_list
%                 face_list
%                 ismember(edge_list(i),IBC_edge)
% %                 if ismember(edge_list(i),IBC_edge)
% %                     if_refine = true;
% %                 else
% %                     if_refine = false;
% %                 end
% %                 if_refine
%             end
%         end
%     end
    
    index = find(~face_list);
    if index
        %edge_list
        %face_list
        %if 
        if ismember(edge_list(index-4),IBC_edge)
            if_refine = true;
        else
            if_refine = false;
        end
    else
        if_refine = true;
    end
    
end
% if ele == 9085
%     edge_list
%     face_list
%     current_level
%     if_refine
% end
end