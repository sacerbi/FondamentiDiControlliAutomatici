function Esempio2_8()
  close all
  clear all
  clc
  % Si analizzi il comportamento del sistema centrifuga descritto dall'esempio
  % Si assuma come stato del sistema x = [w]. Si assuma che w sia pure variabile
  % d'uscita del sistema.
  J = 1;            % Momento di inerzia
  k = 0.1;          % costante caratteristica di coppia d'attrito
  u1 = 0.0;         % coppia motrice nulla
  u2 = 0.5;         % coppia motrice
  xi = -4;          % stato iniziale
  xf = 4;           % stato finale
  w0 = 1;           % Condizione iniziale

  dx = 0.01;
  x = xi:dx:xf;

  w1 = [];
  w2 = [];

  for i = 1 : length(x)
    w1(i) = -(k/J)*(x(i)^2)*sign(x(i)) + u1/J;
    w2(i) = -(k/J)*(x(i)^2)*sign(x(i)) + u2/J;
  endfor

  plot(x,w1,x,w2);
  legend('\omega(t), u(t) = 0','\omega(t), u(t) = 0.5');
  grid;

endfunction

