function output = timeAveragePlot(num_start,nstp,ufile,pfile,meshFile,x,y)
global casenum;
global itaucase
%x = [10.5:0.05:50];
%y = ones(1,length(x))*10;

%% around cylinder
% x = 9.5:0.01:10.5;
% r = 0.5;
% y = sqrt(r^2-(x-10).^2)+10;

%y = [10.5:0.05:50];
%x = ones(1,length(x))*10;

i = 0;
ll = nstp/10;
var = zeros(ll,length(x));

for istp = num_start:nstp+num_start-1
    if mod(istp,10)==0
        i = i+1;
        
        ufilename = strcat(ufile,'100-un-',string(istp+1),'.txt')
        pfilename = strcat(pfile,'100-pn-',string(istp+1),'.txt');
        uv = load(ufilename);
        p = load(pfilename);
        va = [uv,p];
        
        var(i,:) = LinePlot(meshFile,va,x,y,'u');
    end
        
end
%var
meanvar = mean(var);
output = reshape(var,[ll,length(x)]);
 figure(3)
rmsu = zeros(length(y),1);
for i = 1:length(y)
    rmsu(i) = sum(var(:,i).^2);
end
plot(y,rmsu,'b-')
hold on

%var = reshape(var,[],1);
filename = strcat('./timeAverage/static','p3','urmsx.txt');
fileID = fopen(filename,'w');
fprintf(fileID,'%6.20f \n',output);
fclose(fileID);



% filename = strcat('./testvtk/cases/solution/timeAverage/',casenum,itaucase,...
%     'p3ref','urms.txt');
% uref = load(filename);
% uref = reshape(uref,[ll,length(x)]);
% rmsu = zeros(length(x),1);
% for i = 1:length(x)
%     rmsu(i) = sum(uref(:,i).^2);
% end
% plot(x,rmsu,'k-')
% filename = strcat('./testvtk/cases/solution/timeAverage/',itaucase,'p2_fine/',...
%     casenum,'ubar.txt');
% ubar = load(filename);
% ll = 94
% ubar = reshape(ubar,[ll,length(x)]);
% figure(1)
% hold on
% plot(x,mean(ubar),'m-')
% hold on
end
