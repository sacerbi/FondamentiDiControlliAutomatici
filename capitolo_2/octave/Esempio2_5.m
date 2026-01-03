function Esempio2_5()
  close all
  clear all
  clc
  % Si analizzi il comportamento del sistema carrello-molla descritto dall'esempio
  % Si assuma come stato del sistema x = [s sdot]. Si assuma invece come uscita
  % l'energia totale ET.
  Fm = 5;           % Forza applicata al carrello
  kt0 = 10;         % coefficiente elastico a riposo
  alpha = 2;        % decadimento della capacità elastica
  t0 = 0;           % istante iniziale
  s0 = 0;           % Condizione iniziale
  s0dot = 0;        % Condizione iniziale
  x0 = [s0 s0dot];  % Vettore C.I.
  M = 10;           % Massa del carrello

  dt = 0.01;
  t = 0:dt:10;

  [T, x] = ode45(@sistemaCarrelloMolla, t, [x0], M, Fm, kt0, alpha, t0);

  y = (1/2) * (kt0 * e.^(-alpha * (T - t0)) .* x(:,1).^2 + M * x(:,2).^2);

  plot(T,x(:,1),";s(t);",T,x(:,2),";sdot(t);",T,y,";ET(t);");
  grid;

endfunction

function xp = sistemaCarrelloMolla(t, x, M, Fm, kt0, alpha, t0)
  % Ponendo
  %   u(t) = Fm
  %   x(t) = [s(t) sdot(t)]
  %   y(t) = ET(t)
  % risulta
  %   xdot(t) = [sdot(t) (u(t)-kt0*e^(-alpha*(t-t0))*s(t))/M]
  %   ET(t)   = 1/2 * (kt0*e^(-alpha*(t-t0))*s(t)^2 + M*sdot(t)^2)
  u = Fm;
  Fe = kt0 * e^(-alpha*(t-t0)) * x(1);

  xp = [x(2) ...
        (u - Fe)/M];
endfunction

