function Esercizio2_6()
  close all
  clear all
  clc
  % Si consideri il sistema dinamico con equazione di stato
  %   xdot(t) = x(t)^2 + x(t) + u(t)
  % Si determinino le proprietà di stabilità degli stati di equilibrio del
  % sistema al variare di u(t) = ubar <= 1/4

  % Grafico il comportamento al variare di x in situazione di stabilità
  x = -2:0.01:2;
  u1 = 1/8;
  u2 = 1/4;

  y1 = x.^2 + x + u1*ones(size(x));
  y2 = x.^2 + x + u2*ones(size(x));

  plot(x,y1,x,y2)
  hold on
  axis([-2 2 -1 2])
  plot([-2 2], [0 0], 'k')
  legend('u(t)=1/8','u(t)=1/4')
  grid

  % Simulo per due condizioni di x(0) per ubar < 1/4
  % nella prima mi pongo prima del primo punto di stabilità
  % nella seconda mi pongo tra il primo e il secondo punto di stabilità
  t = 0:0.001:10;
  x0_1 = -1;
  x0_2 = -0.3;
  [T1, y11] = ode45(@sistema, t, [x0_1], u1);
  [T2, y12] = ode45(@sistema, t, [x0_2], u1);

  figure
  plot(t,y11,t,y12)
  legend('x(0)_1 = -1', 'x(0)_2 = -0.3')
  grid

endfunction

function [xp] = sistema(t, x, u)
  % Implementa il sistema
  xp = x.^2+x+u;
endfunction

