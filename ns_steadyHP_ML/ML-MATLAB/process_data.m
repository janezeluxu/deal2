function process_data()
% fileNameo = "../output/mesh/lidcavity-100-1mesh.txt";
% fileNameVertexo = "../output/mesh/lidcavity-100-1v_to_e_indices.txt";
% fileNamef = "../output/mesh/lidcavity-100-5mesh.txt";
% fileNameVertexf = "../output/mesh/lidcavity-100-5v_to_e_indices.txt";
% outputfilename = 'prediction_Re100';
% getmeshType(fileNameo,...
%     fileNameVertexo,fileNamef,fileNameVertexf,outputfilename);
% 
% uv = load("../output/solution/lidcavityRe100un-1-.txt");
% p = load("../output/solution/lidcavityRe100pn-1-.txt");
% [input_feature,ele_sizeo] = ...
%     load_data(fileNameo,fileNameVertexo,uv,p);
% save('featureRe100.mat','input_feature')
% 
% fileNameo = "../output/mesh/lidcavity-400-1mesh.txt";
% fileNameVertexo = "../output/mesh/lidcavity-400-1v_to_e_indices.txt";
% fileNamef = "../output/mesh/lidcavity-400-5mesh.txt";
% fileNameVertexf = "../output/mesh/lidcavity-400-5v_to_e_indices.txt";
% outputfilename = 'prediction_Re400';
% getmeshType(fileNameo,...
%     fileNameVertexo,fileNamef,fileNameVertexf,outputfilename);
% 
% uv = load("../output/solution/lidcavityRe400un-1-.txt");
% p = load("../output/solution/lidcavityRe400pn-1-.txt");
% [input_feature,ele_sizeo] = ...
%     load_data(fileNameo,fileNameVertexo,uv,p);
% save('featureRe400.mat','input_feature')
% 
% 
% fileNameo = "../output/mesh/lidcavity-700-1mesh.txt";
% fileNameVertexo = "../output/mesh/lidcavity-700-1v_to_e_indices.txt";
% fileNamef = "../output/mesh/lidcavity-700-5mesh.txt";
% fileNameVertexf = "../output/mesh/lidcavity-700-5v_to_e_indices.txt";
% outputfilename = 'prediction_Re700';
% getmeshType(fileNameo,...
%     fileNameVertexo,fileNamef,fileNameVertexf,outputfilename);
% 
% uv = load("../output/solution/lidcavityRe700un-1-.txt");
% p = load("../output/solution/lidcavityRe700pn-1-.txt");
% [input_feature,ele_sizeo] = ...
%     load_data(fileNameo,fileNameVertexo,uv,p);
% save('featureRe700.mat','input_feature')
% 
% fileNameo = "../output/mesh/lidcavity-1000-1mesh.txt";
% fileNameVertexo = "../output/mesh/lidcavity-1000-1v_to_e_indices.txt";
% fileNamef = "../output/mesh/lidcavity-1000-5mesh.txt";
% fileNameVertexf = "../output/mesh/lidcavity-1000-5v_to_e_indices.txt";
% outputfilename = 'prediction_Re1000';
% cellSize = getmeshType(fileNameo,...
%     fileNameVertexo,fileNamef,fileNameVertexf,outputfilename);
% 
% uv = load("../output/solution/lidcavityRe1000un-1-.txt");
% p = load("../output/solution/lidcavityRe1000pn-1-.txt");
% [input_feature,ele_sizeo] = ...
%     load_data(fileNameo,fileNameVertexo,uv,p);
% save('featureRe1000.mat','input_feature')
% 
fileNameo = "../output/mesh/lidcavity-800-1mesh.txt";
fileNameVertexo = "../output/mesh/lidcavity-800-1v_to_e_indices.txt";
fileNamef = "../output/mesh/lidcavity-800-5mesh.txt";
fileNameVertexf = "../output/mesh/lidcavity-800-5v_to_e_indices.txt";
outputfilename = 'prediction_Re800';
cellSize = getmeshType(fileNameo,...
    fileNameVertexo,fileNamef,fileNameVertexf,outputfilename);

uv = load("../output/solution/lidcavityRe800un-1-.txt");
p = load("../output/solution/lidcavityRe800pn-1-.txt");
[input_feature,ele_sizeo] = ...
    load_data(fileNameo,fileNameVertexo,uv,p);
save('featureRe800.mat','input_feature')
% 
% fileNameo = "../output/mesh/boundary-100-1mesh.txt";
% fileNameVertexo = "../output/mesh/boundary-100-1v_to_e_indices.txt";
% fileNamef = "../output/mesh/boundary-100-5mesh.txt";
% fileNameVertexf = "../output/mesh/boundary-100-5v_to_e_indices.txt";
% outputfilename = 'prediction_boundary_Re100';
% cellSize = getmeshType(fileNameo,...
%     fileNameVertexo,fileNamef,fileNameVertexf,outputfilename);
% 
% uv = load("../output/solution/boundaryRe100un-1-.txt");
% p = load("../output/solution/boundaryRe100pn-1-.txt");
% [input_feature,ele_sizeo] = ...
%     load_data(fileNameo,fileNameVertexo,uv,p);
% save('feature_boundary_Re100.mat','input_feature')
% 
% fileNameo = "../output/mesh/boundary-1000-1mesh.txt";
% fileNameVertexo = "../output/mesh/boundary-1000-1v_to_e_indices.txt";
% fileNamef = "../output/mesh/boundary-1000-5mesh.txt";
% fileNameVertexf = "../output/mesh/boundary-1000-5v_to_e_indices.txt";
% outputfilename = 'prediction_boundary_Re1000';
% cellSize = getmeshType(fileNameo,...
%     fileNameVertexo,fileNamef,fileNameVertexf,outputfilename);
% 
% uv = load("../output/solution/boundaryRe1000un-1-.txt");
% p = load("../output/solution/boundaryRe1000pn-1-.txt");
% [input_feature,ele_sizeo] = ...
%     load_data(fileNameo,fileNameVertexo,uv,p);
% save('feature_boundary_Re1000.mat','input_feature')
% 
% 
% fileNameo = "../output/mesh/backstep-100-1mesh.txt";
% fileNameVertexo = "../output/mesh/backstep-100-1v_to_e_indices.txt";
% fileNamef = "../output/mesh/backstep-100-5mesh.txt";
% fileNameVertexf = "../output/mesh/backstep-100-5v_to_e_indices.txt";
% outputfilename = 'prediction_backstep_Re100';
% cellSize = getmeshType(fileNameo,...
%     fileNameVertexo,fileNamef,fileNameVertexf,outputfilename);
% 
% uv = load("../output/solution/backstepRe100un-1-.txt");
% p = load("../output/solution/backstepRe100pn-1-.txt");
% [input_feature,ele_sizeo] = ...
%     load_data(fileNameo,fileNameVertexo,uv,p);
% save('feature_backstep_Re100.mat','input_feature')

% fileNameo = "../output/mesh/cylinder-50-1mesh.txt";
% fileNameVertexo = "../output/mesh/cylinder-50-1v_to_e_indices.txt";
% fileNamef = "../output/mesh/cylinder-50-5mesh.txt";
% fileNameVertexf = "../output/mesh/cylinder-50-5v_to_e_indices.txt";
% outputfilename = 'prediction_cylinder_Re50';
% cellSize = getmeshType(fileNameo,...
%     fileNameVertexo,fileNamef,fileNameVertexf,outputfilename);
% 
% uv = load("../output/solution/cylinderRe50un-1-.txt");
% p = load("../output/solution/cylinderRe50pn-1-.txt");
% [input_feature,ele_sizeo] = ...
%     load_data(fileNameo,fileNameVertexo,uv,p);
% save('feature_cylinder_Re50.mat','input_feature')
% 
% fileNameo = "../output/mesh/Shape1-1mesh.txt";
% fileNameVertexo = "../output/mesh/Shape1-1v_to_e_indices.txt";
% fileNamef = "../output/mesh/Shape1-5mesh.txt";
% fileNameVertexf = "../output/mesh/Shape1-5v_to_e_indices.txt";
% outputfilename = 'prediction_Shape1_Re50';
% cellSize = getmeshType(fileNameo,...
%     fileNameVertexo,fileNamef,fileNameVertexf,outputfilename);
% 
% uv = load("../output/solution/Shape1un-1-.txt");
% p = load("../output/solution/Shape1pn-1-.txt");
% [input_feature,ele_sizeo] = ...
%     load_data(fileNameo,fileNameVertexo,uv,p);
% save('feature_Shape1_Re50.mat','input_feature')

% fileNameo = "../output/mesh/Shape2-1mesh.txt";
% fileNameVertexo = "../output/mesh/Shape2-1v_to_e_indices.txt";
% fileNamef = "../output/mesh/Shape2-5mesh.txt";
% fileNameVertexf = "../output/mesh/Shape2-5v_to_e_indices.txt";
% outputfilename = 'prediction_Shape2_Re50';
% cellSize = getmeshType(fileNameo,...
%     fileNameVertexo,fileNamef,fileNameVertexf,outputfilename);
% 
% uv = load("../output/solution/Shape2un-1-.txt");
% p = load("../output/solution/Shape2pn-1-.txt");
% [input_feature,ele_sizeo] = ...
%     load_data(fileNameo,fileNameVertexo,uv,p);
% save('feature_Shape2_Re50.mat','input_feature')
% 
% fileNameo = "../output/mesh/Shape3-1mesh.txt";
% fileNameVertexo = "../output/mesh/Shape3-1v_to_e_indices.txt";
% fileNamef = "../output/mesh/Shape3-5mesh.txt";
% fileNameVertexf = "../output/mesh/Shape3-5v_to_e_indices.txt";
% outputfilename = 'prediction_Shape3_Re50';
% cellSize = getmeshType(fileNameo,...
%     fileNameVertexo,fileNamef,fileNameVertexf,outputfilename);
% 
% uv = load("../output/solution/Shape3un-1-.txt");
% p = load("../output/solution/Shape3pn-1-.txt");
% [input_feature,ele_sizeo] = ...
%     load_data(fileNameo,fileNameVertexo,uv,p);
% save('feature_Shape3_Re50.mat','input_feature')

%% load vtk vile and add cellSize at end
%vtkfile = '../output/paraview/Shape1-01.vtk';
vtkfile = '../output/paraview/lidcavityRe800-01.vtk';
fid = fopen(vtkfile, 'at');
if fid ~= -1
    fprintf ( fid, '\n CELL_DATA %d\n', ele_sizeo );
    fprintf ( fid, '\n SCALARS inputcellSize3 int\n' );
    fprintf ( fid, 'LOOKUP_TABLE default\n' );
    
    for element = 1 : ele_sizeo
        fprintf ( fid, '  %d', cellSize(element) );
    end
    fclose(fid);
end
end