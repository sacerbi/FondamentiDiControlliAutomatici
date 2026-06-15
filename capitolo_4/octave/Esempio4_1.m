function Esempio4_1()
  close all
  clear all
  clc
  %
  % Il seguente file produce i grafici presentati nel libro di testo
  %
  alfa1=3;
  alfa2=5;
  punti=-2:0.2:2;
  [X,Y] = meshgrid(punti,punti);
  Z1 = alfa1*X.^2 + alfa2*Y.^2;
  figure;
  surface(X,Y,Z1);
  axis([-2 2 -2 2 0 35]);
  grid;
  xlabel('x_1');
  ylabel('x_2');
  zlabel('W');
  view(3);

  alfa1=1;
  alfa2=1;
  beta=-1;
  Z = alfa1*X.^2 + alfa2*Y.^2+beta*X.^2.*Y;
  figure;
  surface(X,Y,Z);
  axis([-2 2 -2 2 0 10]);
  grid;
  xlabel('x_1');
  ylabel('x_2');
  zlabel('W');
  view(3);
  alfa=1;
  Z = alfa*X.^2;
  figure;
  surface(X,Y,Z);
  axis([-2 2 -2 2 0 4]);
  grid;
  xlabel('x_1');
  ylabel('x_2');
  zlabel('W');
  view(3);
end
