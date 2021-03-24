function [mdl] = ML_training()
input_feature = load("featureRe100.mat");
output1 = load("prediction_Re100");
input1 = input_feature.input_feature;
input_feature = load("featureRe400.mat");
output2 = load("prediction_Re400");
input2 = input_feature.input_feature;
input_feature = load("featureRe700.mat");
output3 = load("prediction_Re700");
input3 = input_feature.input_feature;
input_feature = load("featureRe1000.mat");
output4 = load("prediction_Re1000");
input4 = input_feature.input_feature;

%input_feature = load("featureRe800.mat");
%output5 = load("prediction_Re800");
%input5 = input_feature.input_feature;

% input_feature = load("feature_boundary_Re100.mat");
% output6 = load("prediction_boundary_Re100");
% input6 = input_feature.input_feature;
% input_feature = load("feature_boundary_Re1000.mat");
% output7 = load("prediction_boundary_Re1000");
% input7 = input_feature.input_feature;
% 
% input_feature = load("feature_backstep_Re100.mat");
% output8 = load("prediction_backstep_Re100");
% input8 = input_feature.input_feature;

input_feature = load("feature_cylinder_Re50.mat");
output9 = load("prediction_cylinder_Re50");
input9 = input_feature.input_feature;
%sinput_feature = load("feature_Shape1_Re50.mat");
%output10 = load("prediction_Shape1_Re50");
%input10 = input_feature.input_feature;
input_feature = load("feature_Shape2_Re50.mat");
output11 = load("prediction_Shape2_Re50");
input11 = input_feature.input_feature;
input_feature = load("feature_Shape3_Re50.mat");
output12 = load("prediction_Shape3_Re50");
input12 = input_feature.input_feature;

input = [input1;input2;input3;input4;input9;input11;input12];
output = [output1;output2;output3;output4;output9;output11;output12];

%input = [input9;input11;input12];
%output = [output9;output11;output12];
costMatrix= [0,1,1,1,1;5,0,1,1,1;10,5,0,1,1;20,15,10,0,1;40,30,20,10,0];
% mdl = fitcensemble(input,output,'learners','tree','method','RUSBoost',...
%     'Cost',costMatrix);
mdl = fitcensemble(input,output,'OptimizeHyperparameters','auto',...
    'Cost',costMatrix);
save('trainedModel.mat','mdl');
mdl.resubLoss
confusionchart(output,mdl.resubPredict)
prediction = mdl.resubPredict;


% vertex_size =625;
% fileName = "../output/mesh/lidcavity-1000-1mesh.txt";
% fileNameVertex = "../output/mesh/lidcavity-1000-1v_to_e_indices.txt";
% uv = load("../output/solution/lidcavityRe1000un-1-.txt");
% p = load("../output/solution/lidcavityRe1000pn-1-.txt");
% cellType = load("eleType_response2");
% [input_feature,output,ele_sizeo] = ...
%     load_data(vertex_size,fileName,fileNameVertex,uv,p,cellType);
% prediction = mdl.predict(input_feature);
% confusionchart(output,prediction);
% 
% %% load vtk vile and add cellSize at end
% vtkfile = '../output/paraview/lidcavityRe1000-01.vtk';
% fid = fopen(vtkfile, 'at');
% if fid ~= -1
%     %fprintf(fid, ...... whatever.....
%     %fprintf ( fid, '\n CELL_DATA %d\n', ele_sizeo );  
%     fprintf ( fid, '\n SCALARS responceSize int\n' );
%     fprintf ( fid, 'LOOKUP_TABLE default\n' );
%     
%     for element = 1 : ele_sizeo
%         fprintf ( fid, '  %d', prediction(element) );
%     end
%     fclose(fid);
% end
end