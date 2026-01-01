function Esempio2_3()
  close all
  clear all
  clc
  % Si implementi il circuito dell'esempio e se ne rappresenti l'evoluzione
  R = 1;
  C = 2;
  Vg = 5;
  dt = 0.001;
  t = 0:dt:10;
  x0 = 1;

  [t, y] = ode45(@circuitoRC, t, [x0], R, C, Vg);

  Y = Vg - y;
  plot(t, Y, ";Vr;", t, y, ";Vc;");

endfunction

function [xp] = circuitoRC(t, x, R, C, U)
  % Questa funzione implementa il circuito RC dell'esempio 2.2
  xp = (1/(R*C))*(U - x);
endfunction
