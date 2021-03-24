function [cellSize] = getmeshType(fileNameo,...
    fileNameVertexo,fileNamef,fileNameVertexf,outputfilename)
%vertex_sizeo =625;
%ele_sizeo = 576;
%fileNameo = "../output/mesh/lidcavity-1mesh.txt";
%fileNameVertexo = "../output/mesh/lidcavity-1v_to_e_indices.txt";
[solution_ien_mesho,p_ieno,cellvertexDataxo,cellvertexDatayo,...
    vertexData_mesho,v_e_mesho] = readin_mesh(fileNameo,fileNameVertexo);
%v_e_mesho
ele_sizeo = size(solution_ien_mesho,1);
vertex_sizeo = size(v_e_mesho,1);
%vertex_size = 10442;
%fileNamef = "../output/mesh/lidcavity-5mesh.txt";
%fileNameVertexf = "../output/mesh/lidcavity-5v_to_e_indices.txt";
[solution_ien_meshf,p_ienf,cellvertexDataxf,cellvertexDatayf,...
    vertexData_meshf,v_e_meshf] = readin_mesh(fileNamef,fileNameVertexf);
%v_e_meshf{1}

cellType = cell(ele_sizeo,1);
vertex_sizeo
for i = 1:vertex_sizeo
    v_cood_o = vertexData_mesho(i,:);
    v_cood_f = vertexData_meshf(i,:);
    
    sameNode = testifSameNode(v_cood_o,v_cood_f);
    
    if (sameNode)
        
        eleListo = v_e_mesho{i,1};
        eleListf = v_e_meshf{i,1};
        cell_coord_ox = zeros(length(eleListo),4);
        cell_coord_oy = zeros(length(eleListo),4);
        cell_coord_fx = zeros(length(eleListo),4);
        cell_coord_fy = zeros(length(eleListo),4);
        for j = 1:length(eleListo)
            eleo = eleListo(j);
            cell_coord_ox(j,:) = cellvertexDataxo(eleo,:);
            cell_coord_oy(j,:) = cellvertexDatayo(eleo,:);
            elef = eleListf(j);
            cell_coord_fx(j,:) = cellvertexDataxf(elef,:);
            cell_coord_fy(j,:) = cellvertexDatayf(elef,:);
        end  
        ifsubeleList = zeros(length(eleListo),length(eleListf));
        for j = 1:length(eleListo)
            for k = 1:length(eleListf)
                [ifsubele,in] = testIfsubele(cell_coord_ox(j,:),cell_coord_oy(j,:),...
                    cell_coord_fx(k,:),cell_coord_fy(k,:));
                if ifsubele ==1
                    cell_area_o = getcellarea(cell_coord_ox(j,:),cell_coord_oy(j,:));
                    cell_area_f = getcellarea(cell_coord_fx(k,:),cell_coord_fy(k,:));
                    cell_level = cell_area_o/cell_area_f;
                    cellType{eleListo(j),1} = [cellType{eleListo(j),1},cell_level];
                end
                ifsubeleList(j,k) = ifsubele;
            end
        end
        sumifsub = sum(ifsubeleList,2);
        for j = 1:length(eleListo)
            if sumifsub(j)~=1
                cell_coord_ox
                cell_coord_oy
                cell_coord_fx
                cell_coord_fy
                ifsubeleList
            end
        end
      
                     
        
    else
        disp("not same node, try to find the node")
    end
end

%cellType
cellType
cellSize = zeros(ele_sizeo,1);
for i = 1:ele_sizeo
    cellSize(i) = min(cellType{i,1});
    eleLevel = unique(round(cellType{i,1}));
    if length(eleLevel)>1
        v = round(cellType{i,1});
        [occur,value] = groupcounts(v');
        indices = find(occur == max(occur(:)));
        len = length(indices);
        val_all = zeros(len,1);
        for j = 1:len
            index = indices(j);
            val_all(j) = value(index);
        end
        max_val = max(val_all);
        min_val = min(val_all);
        
        cellSize(i) = max_val;
    end
    
    cellSize(i) = max(cellType{i,1});
end
%cellSize'
cellSize = round(log2(cellSize))/2;
cellSize = floor(cellSize);
cellSize'

fileID = fopen(outputfilename,'w');
formatSpec = '%d \n';
[nrows,ncols] = size(cellSize);
value = zeros(nrows,1);
value(:,1) = cellSize;
for row = 1:nrows
    fprintf(fileID,formatSpec,value(row,:));
end
fclose(fileID);

total_cell_num = 0;
for i = 1:ele_sizeo
    n = cellSize(i);
    total_cell_num = total_cell_num+2^(2*n);
end
total_cell_num
%% load vtk vile and add cellSize at end
% vtkfile = '../output/paraview/lidcavityRe400-01.vtk';
% fid = fopen(vtkfile, 'at');
% if fid ~= -1
%     %fprintf(fid, ...... whatever.....
%     %fprintf ( fid, '\n CELL_DATA %d\n', ele_sizeo );  
%     fprintf ( fid, '\n SCALARS inputcellSize int\n' );
%     fprintf ( fid, 'LOOKUP_TABLE default\n' );
%     
%     for element = 1 : ele_sizeo
%         fprintf ( fid, '  %d', cellSize(element) );
%     end
%     fclose(fid);
% end
end

function cell_area = getcellarea(cell_coord_ox,cell_coord_oy)
x1 = cell_coord_ox(1);
x2 = cell_coord_ox(2);
x4 = cell_coord_ox(3);
x3 = cell_coord_ox(4);

y1 = cell_coord_oy(1);
y2 = cell_coord_oy(2);
y4 = cell_coord_oy(3);
y3 = cell_coord_oy(4);
hh = x1*y2-y1*x2+x2*y3-y2*x3+x3*y4-y3*x4+x4*y1-y4*x1;
cell_area =hh/2;
end

function sameNode = testifSameNode(v_cood_o,v_cood_f)
if (abs(v_cood_o(1)-v_cood_f(1))<1e-6 && abs(v_cood_o(2)-v_cood_f(2))<1e-6)
    sameNode = true;
else
    sameNode = false;
end
 
end

function [ifsubele,in] = testIfsubele(cell_coord_ox_r,cell_coord_oy_r,...
                    cell_coord_fx,cell_coord_fy)
 cell_coord_ox =  [cell_coord_ox_r(1),cell_coord_ox_r(2),cell_coord_ox_r(4),cell_coord_ox_r(3)];
 cell_coord_oy = [cell_coord_oy_r(1),cell_coord_oy_r(2),cell_coord_oy_r(4),cell_coord_oy_r(3)];
 in = zeros(length(cell_coord_ox),1);
 for i = 1:length(cell_coord_ox)  
    in(i) = inpolygon(cell_coord_fx(i),cell_coord_fy(i),cell_coord_ox,cell_coord_oy);
 end
 sumin = sum(in);
 if sumin>2
     ifsubele = 1;
 elseif sumin ==2 %%2 points in element
     %in
     k = find(in);
     index1 = k(1);
     index2 = k(2);
     x_cord = cell_coord_fx(index1);
     y_cord = cell_coord_fy(index1);
     %% find distance to 4 edges
     [dis1] = findDis(x_cord,y_cord,cell_coord_ox,cell_coord_oy);
     x_cord = cell_coord_fx(index2);
     y_cord = cell_coord_fy(index2);
     %% find distance to 4 edges
     [dis2] = findDis(x_cord,y_cord,cell_coord_ox,cell_coord_oy);
     if (min(dis1)>1e-3 || min(dis2)>1e-3)
            ifsubele = 1;
     else
            ifsubele = 0;
     end
        
 elseif sumin ==1 %% find 1 point in element
     %in
     index = find(in);
     %% see if this point is within element
        x_cord = cell_coord_fx(index);
        y_cord = cell_coord_fy(index);
        %% find distance to 4 edges
        [dis] = findDis(x_cord,y_cord,cell_coord_ox,cell_coord_oy);
        if min(dis)>1e-3
            ifsubele = 1;
        else
            ifsubele = 0;
        end
 else
     ifsubele = 0;
 end
end
function [dis] = findDis(x_cord,y_cord,cell_coord_ox,cell_coord_oy)
pt = [x_cord,y_cord,0];
%dis = zeros(length(cell_coord_ox));
coord_x = [cell_coord_ox(1),cell_coord_ox(2),cell_coord_ox(3),cell_coord_ox(4)];
coord_y = [cell_coord_oy(1),cell_coord_oy(2),cell_coord_oy(3),cell_coord_oy(4)];

v1 = [coord_x(1),coord_y(1),0];
v2 = [coord_x(2),coord_y(2),0];
v3 = [coord_x(3),coord_y(3),0];
v4 = [coord_x(4),coord_y(4),0];

dis(1) = point_to_lineseg(pt, v1, v2);
dis(2) = point_to_lineseg(pt, v2, v3);
dis(3) = point_to_lineseg(pt, v3, v4);
dis(4) = point_to_lineseg(pt, v4, v1);
end

function d = point_to_lineseg(pt, v1, v2)
      a = v1 - v2;
      b = pt - v2;
      dist = norm(cross(a,b)) / norm(a);
      
      d_end1 = sqrt((pt(1)-v1(1))^2+(pt(2)-v1(2))^2);
      d_end2 = sqrt((pt(1)-v2(1))^2+(pt(2)-v2(2))^2);
      Line_length = sqrt((v1(1)-v2(1))^2+(v1(2)-v2(2))^2);
      if (max(d_end1,d_end2)<Line_length)
          d = dist;
      else
          d = min(d_end1,d_end2);
      end
          
      
end

function  [ordered_eleo,ordered_ox,ordered_oy] = reorder(cell_coord_ox,cell_coord_oy,eleListo)
min_ox = min(cell_coord_ox,[],2);
[sort_ox, Index_ox] = sort(min_ox);
l = size(min_ox);
if l==1
    ordered_eleo = eleListo;
    ordered_ox = cell_coord_ox;
    ordered_oy = cell_coord_oy;
    
elseif l==2
    cell_coord_ox
    cell_coord_oy
    eleListo

elseif l==4
    xx = zeros(l,1);
    for i = 1:2
        xx(Index_ox(i)) = 0;
    end
    for i = 3:4
        xx(Index_ox(i)) = 1;
    end
    
    min_oy = min(cell_coord_oy,[],2);
    [sort_oy, Index_oy] = sort(min_oy);
    yy = zeros(4,1);
    for i = 1:2
        yy(Index_oy(i)) = 0;
    end
    for i = 3:4
        yy(Index_oy(i)) = 1;
    end
    
    reorder = zeros(4,1);
    for i = 1:4
        if (xx(i) == 0 && yy(i) ==0)
            reorder(i) = 1;
        elseif (xx(i) == 1 && yy(i) == 0)
            reorder(i) = 2;
        elseif (xx(i) == 0 && yy(i) == 1)
            reorder(i) = 3;
        elseif (xx(i) == 1 && yy(i) == 1)
            reorder(i) = 4;
        end
    end
    ordered_eleo = zeros(4,1);
    ordered_ox = zeros(4,4);
    ordered_oy = zeros(4,4);
    for i = 1:4
        ordered_eleo(i) = eleListo(reorder(i));
        ordered_ox(i,:) = cell_coord_ox(reorder(i),:);
        ordered_oy(i,:) = cell_coord_oy(reorder(i),:);
        
    end
else
    disp("y")
end

end