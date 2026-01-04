function Esempio2_9()
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

  dt = 0.001;
  t = 0:dt:1/alpha;

  [T, x] = ode45(@sistemaCarrelloMolla, t, [x0], M, Fm, kt0, alpha, t0);

  y = (1/2) * (kt0 * e.^(-alpha * (T - t0)) .* x(:,1).^2 + M * x(:,2).^2);

  [T1, x1] = ode45(@sistemaCarrelloMollaSemplificato, t, [x0], M, Fm, kt0, alpha, t0);

  y1 = (1/2) * (0.95*kt0*x1(:,1).^2 + M*x1(:,2).^2);

  figure(1)
  plot(T,x(:,1),";s(t);",T,x1(:,1),";s_{teorico}(t);");
  grid;

  figure(2)
  plot(T,x(:,2),";\dot{s}(t);",T,x1(:,2),";\dot{s}_{teorico}(t);");
  grid;

  figure(3)
  plot(T,y,";E_T(t);",T,y1,";E_T_{teorico}(t);");
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

function xp = sistemaCarrelloMollaSemplificato(t, x, M, Fm, kt0, alpha, t0)
  % Ponendo
  %   u(t) = Fm
  %   x(t) = [s(t) sdot(t)]
  %   y(t) = ET(t)
  % risulta
  %   xdot(t) = [sdot(t) (u(t)-kt0*0.95*s(t))/M]
  %   ET(t)   = 1/2 * (kt0*0.95*s(t)^2 + M*sdot(t)^2)
  u = Fm;
  Fe = kt0 * 0.95 * x(1);

  xp = [x(2) ...
        (u - Fe)/M];
endfunction

