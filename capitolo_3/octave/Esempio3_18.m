function Esempio3_18()
  % Comparazione della soluzione reale con il sistema linearizzato.
  close all
  clear all
  clc
  pkg load control
  % Dati
  M=10;
  L=0.5;
  k=2;
  g=9.8;
  x0=[pi/2, 0];   % posizione e velocità iniziali
  u=0;            % coppia

  dt = 0.1;
  t = 0:dt:10;

  % Simulazione
  [t, x] = ode45(@pendolo, t, [x0], u, M, L, k, g);

  % Sistema linearizzato
  A = [0, 1; 0, -k/(M*L*L)];
  B = [0; 1/(M*L*L)];
  C = [0, 0];
  D = 0;
  ulin = u.*ones(size(t));

  sys = ss(A, B, C, D);
  [yl, tl, xl] = lsim(sys, ulin, t, [x0]);



  figure(1)
  hold on
  plot(t, x(:,1), ";T;", t, xl(:,1), ";T_{lin};");
  plot(t, x(:,2), ";W;", t, xl(:,2), ";W_{lin};");
  title("Posizione & Velocità del pendolo")
  grid on

  y = 1/2 * M*L*L .* x(:,2).^2 - M*g*L.*cos(x(:,1));
  figure(2)
  plot(t,y,";E;", t, yl, ";E_{lin};")
  title("Andamento energia")

endfunction

function [xp] = pendolo(t, x, u, M, L, k, g)
  % Evoluzione del sistema pendolo secondo le equazioni 2.36-2.38
  x1p = x(2);
  x2p = -(g/L)*sin(x(1)) - (k/(M*L*L))*x(2) + 1/(M*L*L)*u;
  xp = [x1p, x2p];
endfunction

