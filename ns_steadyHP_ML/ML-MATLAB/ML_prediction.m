function ML_prediction()
%[mdl] = ML_training();
%save('trainedModel.mat','mdl');
model = load('trainedModel.mat');
mdl = model.mdl;

%input_feature = load("feature_Shape1_Re50.mat");
%output = load("prediction_Shape1_Re50");

input_feature = load("featureRe800.mat");
output = load("prediction_Re800");

prediction = mdl.predict(input_feature.input_feature);
loss = mdl.loss(input_feature.input_feature,output)
confusionchart(output,prediction);

filename = "ML_prediction_lidcavity_Re800";
fileID = fopen(filename,'w');
for ele_count = 1:length(prediction)
    fprintf(fileID,'%d\n',prediction(ele_count));
end

ele_sizeo = size(output,1)
%% load vtk vile and add cellSize at end
vtkfile = '../output/paraview/Shape1-01.vtk';
%vtkfile = '../output/paraview/lidcavityRe800-01.vtk';

fid = fopen(vtkfile, 'at');
if fid ~= -1
    %fprintf(fid, ...... whatever.....
    %fprintf ( fid, '\n CELL_DATA %d\n', ele_sizeo );  
    fprintf ( fid, '\n SCALARS responce4 int\n' );
    fprintf ( fid, 'LOOKUP_TABLE default\n' );
    
    for element = 1 : ele_sizeo
        fprintf ( fid, '  %d', prediction(element) );
    end
    fclose(fid);
end
end