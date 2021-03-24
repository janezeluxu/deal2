function [alg_celllevel,multicycle] = AlgebraicModel()
error = load("../output/mesh/estimated_error_per_cell1.txt");
multicycle = load("prediction_Re800");
Nele = length(error);
count1 = 0;
count2 = 0;
count3 = 0;
count4 = 0;
for i = 1:Nele
    if multicycle(i) == 1
        count1 = count1+1;
    elseif multicycle(i) == 2
         count2 = count2+1; 
    elseif multicycle(i) == 3
         count3 = count3+1; 
    elseif multicycle(i) == 4
         count4 = count4+1; 
     
    end
end

count1
count2
count3
count4
count0 = Nele-count1-count2-count3-count4
[err,ind] = sort(error);

alg_celllevel = ones(Nele,1)*-1;
for i = 1:count0
    alg_celllevel(ind(i)) = 0;
end
for i = count0+1:count0+count1
    alg_celllevel(ind(i)) = 1;
end
for i = count0+count1+1:count0+count1+count2
    alg_celllevel(ind(i)) = 2;
end
for i = count0+count1+count2+1:count0+count1+count2+count3
    alg_celllevel(ind(i)) = 3;
end
for i = count0+count1+count2+count3+1:count0+count1+count2+count3+count4
    alg_celllevel(ind(i)) = 4;
end

filename = "ML_alg_lidcavity_Re800";
fileID = fopen(filename,'w');
for ele_count = 1:length(alg_celllevel)
    fprintf(fileID,'%d\n',alg_celllevel(ele_count));
end

confusionchart(multicycle,alg_celllevel);
alg_celllevel'
ele_sizeo = Nele;
vtkfile = '../output/paraview/lidcavityRe800-01.vtk';
fid = fopen(vtkfile, 'at');
if fid ~= -1
    %fprintf ( fid, '\n CELL_DATA %d\n', ele_sizeo );
    fprintf ( fid, '\n SCALARS alg_celllevel int\n' );
    fprintf ( fid, 'LOOKUP_TABLE default\n' );
    
    for element = 1 : ele_sizeo
        fprintf ( fid, '  %d', alg_celllevel(element) );
    end
    fclose(fid);
end
end