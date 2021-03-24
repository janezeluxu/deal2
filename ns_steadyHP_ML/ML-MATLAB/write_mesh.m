fileName = "../output/mesh/Shape1-1mesh.txt";
fileNameVertex = "../output/mesh/Shape1-1v_to_e_indices.txt";
cellLevel = load("prediction_Shape1_Re50");
IENfilename = "../output/mesh/Shape1-1EdgeInfo.txt";
bcFile = "../output/mesh/Shape1-1bc.txt";

outmeshfile = "../NS-steadyHO-quad/mesh/Shape1.txt";
outbcfile = "../NS-steadyHO-quad/mesh/Shape1bc.txt";
create_mesh(fileName,fileNameVertex,cellLevel,IENfilename,bcFile,outmeshfile,outbcfile)