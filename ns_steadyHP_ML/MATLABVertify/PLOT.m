Re = 100;
figure(1)
[uG,xG,yG] = Gia(Re,0.5,0.0);
plot(yG,uG,'ks','MarkerSize',10)
hold on
y = 0:0.01:1;
[u,v] = Re100();
plot(y,u)

figure(2)
[uG,xG,yG] = Gia(Re,0.0,0.5);
plot(xG,uG,'ks','MarkerSize',10)
hold on
plot(y,v)

figure(3)
Re = 1000;
[uG,xG,yG] = Gia(Re,0.5,0.0);
plot(yG,uG,'ks','MarkerSize',10)
hold on
y = 0:0.01:1;
[u,v] = Re1000();
plot(y,u)

figure(4)
[uG,xG,yG] = Gia(Re,0.0,0.5);
plot(xG,uG,'ks','MarkerSize',10)
hold on
plot(y,v)

figure(5)
Re = 400;
[uG,xG,yG] = Gia(Re,0.5,0.0);
plot(yG,uG,'ks','MarkerSize',10)
hold on
y = 0:0.01:1;
[u,v] = Re400();
plot(y,u)

figure(6)
[uG,xG,yG] = Gia(Re,0.0,0.5);
plot(xG,uG,'ks','MarkerSize',10)
hold on
plot(y,v)