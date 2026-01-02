function Esempio2_4()
  pkg load symbolic
  close all
  clear all
  clc
  % Si implementi il circuito dell'esempio e se ne rappresenti l'evoluzione
  R = 1;
  C = 1;
  Vg = 2;
  w = 3;
  dt = 0.001;
  t = 0:dt:20;
  x0 = 10;
  u = Vg * sin(w*t);

  [t, y] = ode45(@circuitoRC, t, [x0], R, C, Vg, w);

  Y = u' - y;
  plot(t, Y, ";Vr;", t, y, ";Vc;");

  figure;
  gamma = atan(w*R*C);
  xteorico = e.^(-t/(R*C))*x0 + (Vg*w*R*C)/(1+w^2*R^2*C^2)*e.^(-t/(R*C))+((Vg)/(sqrt(1+w^2*R^2*C^2)))*sin(w*t -gamma);
  yteorico = -e.^(-t/(R*C))*x0 - (Vg*w*R*C)/(1+w^2*R^2*C^2)*e.^(-t/(R*C))+((Vg*w*R*C)/(sqrt(1+w^2*R^2*C^2)))*cos(w*t -gamma);
  plot(t, yteorico, ";Vr;", t, xteorico, ";Vc;");

endfunction

function [xp] = circuitoRC(t, x, R, C, Vg, w)
  % Questa funzione implementa il circuito RC dell'esempio 2.2
  % La funzione di input è u(t) = Vg * sin(w*t)
  U = Vg * sin(w*t);
  xp = (1/(R*C))*(U - x);
endfunction
