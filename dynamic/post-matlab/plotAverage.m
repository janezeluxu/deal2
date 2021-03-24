num_start = 3541;
nstp = 560;
%nstp = 560+280;
%nstp = 560+280*2;
ufilename = "../solution/";
pfilename = "../solution/";
meshFile =  "../mesh/p3mesh.txt";
%bcele_file = "../mesh/bcEle.txt";
%y = [0:0.05:20];
%x = ones(1,length(y))*11;
x = [10.5:0.05:50];
y = ones(1,length(x))*10;
%timeAveragePlot(num_start,nstp,ufilename,pfilename,meshFile,x,y)

figure(1)
filename = strcat('./timeAverage/static','p3','urmsx.txt');
uref = load(filename);
ll = length(uref)/length(x);
uref = reshape(uref,[ll,length(x)]);
rmsu = zeros(length(x),1);
for i = 1:length(x)
    rmsu(i) = sum(uref(:,i).^2)/ll;
end
plot(x,rmsu,'b-','LineWidth',2)

% ly = length(y)
% h = plot(rmsu(1:(ly-1)/2+1),y(1:(ly-1)/2+1),'b-','LineWidth',2)
% h.Color =  [165,42,42]/256;
% hold on
% 
% ly = length(y)
% rmsuupper = rmsu((ly-1)/2+1:end);
% size(rmsuupper)
% size(y(1:(ly-1)/2))
% h = plot(rmsu((ly-1)/2+1:end),10-y(1:(ly-1)/2+1),'b-','LineWidth',2)
% h.Color =  [165,42,42]/256;
% hold on

% hold on
% filename = strcat('./timeAverage/static','p3','urmsy2.txt');
% uref = load(filename);
% ll = length(uref)/length(x);
% uref = reshape(uref,[ll,length(x)]);
% rmsu = zeros(length(x),1);
% for i = 1:length(x)
%     rmsu(i) = sum(uref(:,i).^2)/ll;
% end
% %plot(y,rmsu,'r--','LineWidth',2)
% 
% ly = length(y)
% h = plot(rmsu(1:(ly-1)/2+1),y(1:(ly-1)/2+1),'r--','LineWidth',2)
% %h.Color =  [165,42,42]/256;
% hold on
% 
% ly = length(y)
% rmsuupper = rmsu((ly-1)/2+1:end);
% size(rmsuupper)
% size(y(1:(ly-1)/2))
% h = plot(rmsu((ly-1)/2+1:end),10-y(1:(ly-1)/2+1),'r--','LineWidth',2)
% %h.Color =  [165,42,42]/256;
% hold on
% 
% hold on
% filename = strcat('./timeAverage/static','p3','urmsy3.txt');
% uref = load(filename);
% ll = length(uref)/length(x);
% uref = reshape(uref,[ll,length(x)]);
% rmsu = zeros(length(x),1);
% for i = 1:length(x)
%     rmsu(i) = sum(uref(:,i).^2)/ll;
% end
% %plot(y,rmsu,'m-.','LineWidth',2)
% ly = length(y)
% h = plot(rmsu(1:(ly-1)/2+1),y(1:(ly-1)/2+1),'m.','LineWidth',2)
% %h.Color =  [165,42,42]/256;
% hold on
% 
% ly = length(y)
% rmsuupper = rmsu((ly-1)/2+1:end);
% size(rmsuupper)
% size(y(1:(ly-1)/2))
% h = plot(rmsu((ly-1)/2+1:end),10-y(1:(ly-1)/2+1),'m.','LineWidth',2)
% %h.Color =  [165,42,42]/256;
% hold on
% 
% figure(2)
% ly = length(y)
% h = plot(rmsu(1:(ly-1)/2+1),y(1:(ly-1)/2+1),'b-','LineWidth',2)
% h.Color =  [165,42,42]/256;
% hold on
% 
% ly = length(y)
% rmsuupper = rmsu((ly-1)/2+1:end);
% size(rmsuupper)
% size(y(1:(ly-1)/2))
% h = plot(rmsu((ly-1)/2+1:end),10-y(1:(ly-1)/2+1),'r--','LineWidth',2)
% h.Color =  [165,42,42]/256;
% hold on