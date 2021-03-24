function test()
cell_coord_ox =[

  -22.7661  -24.2180  -22.7454  -24.2057
  -22.7454  -24.2057  -22.7368  -24.2021
  -25.6645  -25.6564  -24.2180  -24.2057
  -25.6564  -25.6548  -24.2057  -24.2021];


cell_coord_oy =[

  -19.7068  -19.6771  -21.1612  -21.1394
  -21.1612  -21.1394  -22.6282  -22.6110
  -19.6354  -21.1070  -19.6771  -21.1394
  -21.1070  -22.5853  -21.1394  -22.6110];

cell_coord_ox =[

   -0.1913         0   -0.2317   -0.0329
    0.1913    0.1744         0   -0.0329];


cell_coord_oy =[

    0.4619    0.5000    0.5842    0.6296
    0.4619    0.6044    0.5000    0.6296];

cell_coord_ox =[

    0.4511   -1.3348    0.6325   -1.1160
    0.6325   -1.1160    0.8140   -0.8973
    0.6325    0.8140    2.4132    2.5421
    0.4511    0.6325    2.2843    2.4132];


cell_coord_oy =[

  -19.5848  -20.0362  -21.4679  -21.7521
  -21.4679  -21.7521  -23.3509  -23.4679
  -21.4679  -23.3509  -21.5732  -23.2964
  -19.5848  -21.4679  -19.8500  -21.5732];


cell_coord_fx =[

    0.5418   -0.3418    0.6325   -0.2417
    0.6325   -0.2417    0.7233   -0.1417
    0.6325    0.7233    1.5229    1.6005
    0.5418    0.6325    1.4453    1.5229];


cell_coord_fy =[

  -20.5263  -20.7102  -21.4679  -21.6100
  -21.4679  -21.6100  -22.4094  -22.5097
  -21.4679  -22.4094  -21.5205  -22.4221
  -20.5263  -21.4679  -20.6189  -21.5205];
    
plot(cell_coord_ox(1,:),cell_coord_oy(1,:),'b*')
hold on
plot(cell_coord_ox(2,:),cell_coord_oy(2,:),'r*')
hold on
plot(cell_coord_ox(3,:),cell_coord_oy(3,:),'k*')
hold on
plot(cell_coord_ox(4,:),cell_coord_oy(4,:),'g*')
hold on

plot(cell_coord_fx(1,:),cell_coord_fy(1,:),'b*')
hold on
plot(cell_coord_fx(2,:),cell_coord_fy(2,:),'r*')
hold on
plot(cell_coord_fx(3,:),cell_coord_fy(3,:),'k*')
hold on
plot(cell_coord_fx(4,:),cell_coord_fy(4,:),'g*')
hold on

[ifsubele,in] = testIfsubele(cell_coord_ox(4,:),cell_coord_oy(4,:),...
            cell_coord_fx(3,:),cell_coord_fy(3,:))
        
for j = 1:4
    for k = 1:4
        j
        k
        [ifsubele,in] = testIfsubele(cell_coord_ox(j,:),cell_coord_oy(j,:),...
            cell_coord_fx(k,:),cell_coord_fy(k,:));
        ifsubele
        ifsubeleList(j,k) = ifsubele;
    end
end
ifsubeleList              
%plot(cell_coord_fx,cell_coord_fy,'ro')
%reorder(cell_coord_ox,cell_coord_oy)
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

function  reorder(cell_coord_ox,cell_coord_oy)
min_ox = min(cell_coord_ox,[],2);
[sort_ox, Index_ox] = sort(min_ox);
xx = zeros(4,1);
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
reorder

end