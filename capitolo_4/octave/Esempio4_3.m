function Esempio4_1()
  close all
  clear all
  clc
  %
  % Il seguente file produce i grafici presentati nel libro di testo
  %
  alfa=1;
  beta=1;
  punti=-4:0.2:4;
  [X,Y] = meshgrid(punti,punti);
  Z1 = alfa*(1 - cos(X)) + beta*Y.^2;
  figure;
  surface(X,Y,Z1);
  axis([-4 4 -4 4 0 6]);
  grid;
  xlabel('x_1');
  ylabel('x_2');
  zlabel('W');
  view(3);

  punti=-1:0.1:1;
  [X,Y] = meshgrid(punti,punti);
  Z2 = exp(X.^2 + Y.^2) - 1;
  figure;
  surface(X,Y,Z2);
  axis([-1 1 -1 1 0 7]);
  grid;
  xlabel('x_1');
  ylabel('x_2');
  zlabel('W');
  view(3);

  Z3 = abs(X) + abs(Y);
  figure;
  surface(X,Y,Z3);
  axis([-1 1 -1 1 0 2]);
  grid;
  xlabel('x_1');
  ylabel('x_2');
  zlabel('W');
  view(3);

  Z4 = max(abs(X), abs(Y));
  figure;
  surface(X,Y,Z4);
  axis([-1 1 -1 1 0 1]);
  grid;
  xlabel('x_1');
  ylabel('x_2');
  zlabel('W');
  view(3);
end
