function Bezier_curve()
%Hello! this will plot Bezier curve for n control points 
%This is a replacement of the program 'Parametic Cubic Bezier Curve'
%submitted before ...Bezier for any number of points ...enjoy 
clear all 
clc

% n=input('Enter no. of points  ');
% w=input('Press 1 for entry through mouse or 2 for keyboard repectively-->');
% if w==1
%     axis([-0.5 0.5 -0.5 0.5])
%     [p]=ginput(n);
% end
% if w==2
%     [p]=input('Enter co-odinates of points within brackets ->[x1 y1;x2 y2;x3 y3;...;xn yn] ');
% end
% n
% p

% n = 11;
% p = [    0.4902    0.0029
%     0.3704    0.2569
%     0.1630    0.4029
%     0.0294    0.4905
%    -0.2033    0.4117
%    -0.3715    0.1489
%    -0.3301   -0.1255
%    -0.1388   -0.3299
%    -0.0075   -0.3328
%     0.2598   -0.2774
%     0.4902    0.0029];


n =17;


p =[

    0.49    0.0000
    0.4395    0.1635
    0.3105    0.2453
    0.1630    0.2832
    0.0040    0.3007
   -0.1204    0.2774
   -0.2540    0.2277
   -0.3531    0.1255
   -0.4706    0.0029
   -0.4061   -0.0934
   -0.3070   -0.1839
   -0.1826   -0.2365
    0.0040   -0.3007
    0.1285   -0.2599
    0.2690   -0.1810
    0.3957   -0.0905
    0.49   -0.0];
    
n =17;


% p =[
% 
%     0.20    0.00
%     0.1907    0.1810
%     0.1354    0.3007
%     0.0778    0.3883
%     0.0040    0.4876
%    -0.1250    0.3387
%    -0.1941    0.1956
%    -0.1526    0.0934
%    -0.1987    0.0029
%    -0.2056   -0.1255
%    -0.1526   -0.2190
%    -0.0997   -0.3007
%    -0.0006   -0.4672
%     0.0478   -0.3854
%     0.1285   -0.2920
%     0.1653   -0.1460
%     0.2    0.0000];
 BezierCurve(n,p)
end
 function   []= BezierCurve(n,p)
 n1=n-1;
for    i=0:1:n1
sigma(i+1)=factorial(n1)/(factorial(i)*factorial(n1-i));  % for calculating (x!/(y!(x-y)!)) values 
end
l=[];
UB=[];
for u=0:0.002:1
for d=1:n
UB(d)=sigma(d)*((1-u)^(n-d))*(u^(d-1));
end
l=cat(1,l,UB);                                      %catenation 
end
P=l*p;
l = line(P(:,1),P(:,2))
xx = l.XData';
yy = l.YData';
save('x.mat','xx');
save('y.mat','yy');
%lp = line(p(:,1),p(:,2))

 end