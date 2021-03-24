function plotForce() 

% nstp = 160;
% num_start = 3521;
% fileName =  "../mesh/p5mesh.txt";
% bcele_file = "../mesh/bcEle.txt";
% 
% ufilename = "../solutioncylinderp5/";
% pfilename = "../solutioncylinderp5/";
% p = 5;

% nstp = 10;
% num_start = 3201;
% fileName =  "../mesh/p3mesh.txt";
% bcele_file = "../mesh/bcEle.txt";
% 
% ufilename = "../solutioncylinderp3/";
% pfilename = "../solutioncylinderp3/";
% p = 3;

nstp = 1580;
num_start = 3201;
fileName =  "../mesh/p3mesh.txt";
bcele_file = "../mesh/bcEle.txt";

ufilename = "../solution/";
pfilename = "../solution/";
p = 3;

[Force]=force(p,num_start,nstp,fileName,bcele_file,ufilename,pfilename);
Force(2:end,:)

figure(1)
%plot(num_start:10:nstp+num_start-1,2*ForceDynamic(1:end,1),'r')
%hold on
plot(num_start+10:10:nstp+num_start-1,2*Force(2:end,1),'r')
%ylim([1.2 1.4])

figure(2)
%plot(num_start:10:nstp+num_start-1,2*ForceDynamic(1:end,2),'r')
%hold on
plot(num_start+10:10:nstp+num_start-1,2*Force(2:end,2),'r')
%ylim([-0.5 0.5])
end
