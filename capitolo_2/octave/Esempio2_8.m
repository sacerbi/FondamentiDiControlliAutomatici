function Esempio2_8()
  close all
  clear all
  clc
  % Si analizzi il comportamento del sistema centrifuga descritto dall'esempio
  % Si assuma come stato del sistema x = [w]. Si assuma che w sia pure variabile
  % d'uscita del sistema.
  J = 1;            % Momento di inerzia
  k = 0.1;          % costante caratteristica di coppia d'attrito
  u = 0.5;          % coppia motrice
  t0 = 0;           % istante iniziale
  w0 = 10;           % Condizione iniziale

  dt = 0.01;
  t = t0:dt:10;

  [T, x] = ode45(@sistemaCentrifuga, t, [w0], J, k, u);

  y = x;

  plot(T,y);
  legend('\omega(t)');
  grid;

  % A coppia nulla u(t) = 0
  [T1, x1] = ode45(@sistemaCentrifuga, t, [w0], J, k, 0);

  y1 = x1;
  y1teorico = (J*w0)./(J + k*w0*sign(w0).*t);
  figure
  plot(T1,y1,t,y1teorico);
  legend('\omega(t)', '\omega(t)_{teorico}')
  grid;

  error = (y1 - y1teorico')'*(y1 - y1teorico')

endfunction

function xp = sistemaCentrifuga(t, x, J, k, u)
  % Ponendo
  %   u(t) nota
  %   x(t) = [w(t)]
  %   y(t) = w(t)
  % risulta
  %   xdot(t) = -(k/J)*x(t)^2*sgn(x(t)) + (1/J)*u(t)
  %   y(t)    = x(t)
  segno = sign(x);

  xp = -(k/J)*x*x*segno + (1/J)*u;

endfunction

