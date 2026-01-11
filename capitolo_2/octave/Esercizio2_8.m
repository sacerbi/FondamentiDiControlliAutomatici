function Esercizio2_8()
  close all
  clear all
  clc
  % Considerando il sistema dinamico
  %   xdot(t) = |x(t)| + u(t)
  % Si calcoli il movimento dello stato per ubar=0 e x(0)=x0.
  % Poi per u(t) = ubar si determinino gli stati di equilibrio e le corrispondenti
  % proprietà di stabilità.
  dt = 0.001;
  t0 = 0; tf = 3;
  t = t0:dt:tf;
  x0_1 = 1;
  x0_2 = -1;
  ubar = 0;
  [t1,x1] = ode45(@sistema, t, [x0_1], ubar);
  [t2,x2] = ode45(@sistema, t, [x0_2], ubar);
  figure
  plot(t,x1,t,x2,[t0 tf], [0 0], 'k')
  xlabel('t')
  ylabel('x(t)')
  grid
  title('Esercizio 2.8 - Movimento stato')
  legend('u=0, x(0)=1','u=0, x(0)=-1')

  % Si analizzano ora gli stati di equilibrio ponedo 0 = |xbar| + ubar
  %   - ubar > 0 --> No stati equilibrio
  %   - ubar = 0 --> xbar = 0 equilibrio instabile
  %   - ubar < 0 --> xbar = -ubar instabile e xbar = ubar asintoticamente stabile
  x = -3:0.01:3;
  ubar1 = 1;
  ubar2 = 0;
  ubar3 = -1;
  xdot1 = abs(x) + ubar1*ones(size(x));
  xdot2 = abs(x) + ubar2*ones(size(x));
  xdot3 = abs(x) + ubar3*ones(size(x));
  figure
  plot(x,xdot1,x,xdot2,x,xdot3,[-3 3], [0 0], 'k')
  xlabel('x')
  ylabel('dx/dt')
  grid
  title('Esercizio 2.8 - Studio equilibrio')
  legend('u=1','u=0','u=-1')

  tf = 20;
  t = t0:dt:tf;
  x0_3 = -0.1;
  [t3, x3] = ode45(@sistema, t, [x0_3], ubar3);
  figure
  plot(t,x3,[t0 tf], [0 0], 'k')
  xlabel('t')
  ylabel('x(t)')
  grid
  title('Esercizio 2.8 - Studio equilibrio')
  legend('u=-1, x(0)=-0.1')
endfunction

function [xp] = sistema(t,x,u)
  xp = [abs(x) + u];
endfunction

