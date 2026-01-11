function Esercizio2_7()
  close all
  clear all
  clc
  % Si determinino gli stati di equilibrio e le corrispondneti proprietà di
  % stabilità per i sistemi di ordine 1 con le equazioni di stato seguenti al
  % variare di u(t) = ubar

  % Primo punto
  Esercizio2_7_1()

  % Secondo punto
  Esercizio2_7_2()

  % Terzo punto
  Esercizio2_7_3()

endfunction

function Esercizio2_7_1()
  clear all
  % Per determinare gli stati di equilibrio, analizzo al variare di x
  % Poneno 0 = ubar * xbar^3
  %   - per u = 0 tutti gli stati sono in equilibrio stabile
  %   - per u < 0, x = 0 è punto di equilibrio. Dato il segno di xdot sappiamo
  %     che è globalmente asintoticamente stabile
  %   - per u > 0, x = 0 è punto di equilibrio instabile
  dx = 0.01;
  x = -3:dx:3;
  u1 = -1*ones(size(x));
  u2 = zeros(size(x));
  u3 = 1*ones(size(x));
  xdot1 = u1.*(x.^3);
  xdot2 = u2.*(x.^3);
  xdot3 = u3.*(x.^3);
  figure
  plot(x,xdot1,';u=-1;',x,xdot2,';u=0;',x,xdot3,';u=1;');
  grid
  xlabel('x')
  ylabel('dx/dt')
  title('Esercizio 2.7.1')

endfunction

function Esercizio2_7_2()
  clear all
  % Per determinare gli stati di equilibrio, analizzo al variare di x
  dx = 0.01;
  x_neg = -3:dx:0;
  x_pos = dx:dx:3;
  x = [x_neg x_pos];
  u1_neg = -1*ones(size(x_neg));
  u2_neg = zeros(size(x_neg));
  u3_neg = 1*ones(size(x_neg));
  u4_neg = 2.1*ones(size(x_neg));
  u1_pos = -1*ones(size(x_pos));
  u2_pos = zeros(size(x_pos));
  u3_pos = 1*ones(size(x_pos));
  u4_pos = 2.1*ones(size(x_pos));

  xdot1 = [u1_neg+x_neg, u1_pos-2*(1-exp(-x_pos))];
  xdot2 = [u2_neg+x_neg, u2_pos-2*(1-exp(-x_pos))];
  xdot3 = [u3_neg+x_neg, u3_pos-2*(1-exp(-x_pos))];
  xdot4 = [u4_neg+x_neg, u4_pos-2*(1-exp(-x_pos))];

  figure
  plot(x,xdot1,';u=-1;',x,xdot2,';u=0;',x,xdot3,';u=1;',x,xdot4,';u=2.1;',[-3 3], [0 0], 'k');
  grid
  xlabel('x')
  ylabel('dx/dt')
  title('Esercizio 2.7.2 - Equilibrio')
  % Dalla formulazione di xdot posso intuire che il numero di punti di stabilità
  % varia a seconda del valore di u(t)=ubar.
  % Per u < 0 non ci sono stati di equilibrio.
  % Per u = 0, x=0 è punto di equilibrio instabile perché xdot < 0 nell'intorno
  % Per u > 0, per la formulazione di xdot nei casi di x>0:
  %   - se u > 2, c'è un solo punto di equilibrio x=-u instabile (xdot stesso
  %     segno sia prima che dopo)
  %   - se u in (0,2), ci sono due punti di eq: x=-u instabile e x=-ln(1-1/2*u)
  %     asintoticamente stabile
  % Studio quindi le proprietà di stabilità per u=1, per x(0) < 0 ma maggiore di
  % -u e x(0) > 0
  ubar = 1;
  x0_1 = -0.9;
  x0_2 = 3;
  t = 0:0.001:10;
  xbar = -log(1-ubar/2);
  [t1,xp1] = ode45(@sistema2, t, [x0_1], ubar);
  [t2,xp2] = ode45(@sistema2, t, [x0_2], ubar);
  figure
  plot(t,xp1,';u=1, x(0)=-0.9;',t,xp2,';u=1, x(0)=3;',[0 10], [xbar xbar], 'k');
  grid
  xlabel('t')
  ylabel('x(t)')
  title('Esercizio 2.7.2 - Stabilità')
endfunction

function Esercizio2_7_3()
  clear all
  % Per determinare gli stati di equilibrio, analizzo al variare di x
  % All'equilibrio devo risolvere 0 = xbar(alpha -ubar^2*xbar^2)
  % I casi da analizzare sono:
  %   - ubar = 0 e alpha != 0 --> unico punto in xbar = 0. La stabilità dipende
  %     da alpha
  %   - ubar = 0 e alpha = 0 --> tutti gli stati sono eq stabile
  %   - ubar != 0 e alpha = 0 --> degenera nel sistema 2.7.1
  %   - ubar != 0 e alpha < 0 --> unico punto in xbar = 0 asintoticamente stabile
  %     perché xdot/x < 0
  %   - ubar != 0 e alpha > 0 --> gli stati sono xbar = 0, xbar = +- sqrt(alpha/u^2)
  %     il primo instabile, gli altri asintoticamente stabili
  dx = 0.01;
  x = -3:dx:3;

  % CASO 1
  ubar = 0;
  alpha1 = 1;
  alpha2 = -1;
  xdot1 = alpha1.*x - ((ubar*ones(size(x))).^2).*(x.^3);
  xdot2 = alpha2.*x - ((ubar*ones(size(x))).^2).*(x.^3);
  figure
  plot(x,xdot1,x,xdot2,[-3 3], [0 0], 'k')
  xlabel('x')
  ylabel('dx/dt')
  grid
  legend('u=0, \alpha=1','u=0, \alpha=-1')
  title('Esercizio 2.7.3 - CASO 1')

  % CASO 2
  ubar = 0;
  alpha = 0;
  xdot = alpha.*x - ((ubar*ones(size(x))).^2).*(x.^3);
  figure
  plot(x,xdot,[-3 3], [0 0], 'k')
  xlabel('x')
  ylabel('dx/dt')
  grid
  legend('u=0, \alpha=0')
  title('Esercizio 2.7.3 - CASO 2')

  % CASO 3
  ubar = 1;
  alpha = 0;
  xdot = alpha.*x - ((ubar*ones(size(x))).^2).*(x.^3);
  figure
  plot(x,xdot,[-3 3], [0 0], 'k')
  xlabel('x')
  ylabel('dx/dt')
  grid
  legend('u=1, \alpha=0')
  title('Esercizio 2.7.3 - CASO 3')

  % CASO 4
  ubar = 1;
  alpha = -1;
  xdot = alpha.*x - ((ubar*ones(size(x))).^2).*(x.^3);
  figure
  plot(x,xdot,[-3 3], [0 0], 'k')
  xlabel('x')
  ylabel('dx/dt')
  grid
  legend('u=1, \alpha=-1')
  title('Esercizio 2.7.3 - CASO 4')

  % CASO 5
  ubar = 1;
  alpha = 1;
  xdot = alpha.*x - ((ubar*ones(size(x))).^2).*(x.^3);
  figure
  plot(x,xdot,[-3 3], [0 0], 'k')
  xlabel('x')
  ylabel('dx/dt')
  grid
  legend('u=1, \alpha=1')
  title('Esercizio 2.7.3 - CASO 5')

endfunction

function [xp] = sistema2(t, x, u)
  if x <= 0
    xp = [u + x];
  else
    xp = [u - 2*(1-exp(-x))];
  endif
endfunction

