function Esercizio3_7()
  % Comparazione della soluzione reale con il sistema linearizzato.
  close all
  clear all
  clc
  pkg load control
  % I dati vengono calcolati dalle relazioni di linearizzazione
  % A = df/dx|xbar,ubar
  % B = df/du|xbar,ubar
  % C = dg/dx|xbar,ubar
  % D = dg/du|xbar,ubar
  % Quindi si valuta prima per ubar=1 la condizione di equilibrio:
  % 0 = xbar^2 - ubar*xbar - 2*ubar
  % Si trovano le soluzioni xbar=-1 (ybar = 0) e xbar=2 (ybar=9)
  % Per la seconda soluzione A = 3 --> Non stabile
  A1 = -3;
  B1 = -1;
  C1 = 3;
  D1 = 3;

  ubar = 1;

  dt = 0.1;
  t = 0:dt:10;

  x0 = -1; A = 0.1; w = 2;

  % Sistema linearizzato
  ulin = A*cos(w*t);
  sys = ss(A1, B1, C1, D1);
  yl = lsim(sys, ulin, t);
  % Simulazione sinusoidale
  [teq, xsin] = ode45(@sistema, t, x0, A, w, ubar);
  ysin=xsin.^3+(1+A*cos(w.*t')).^3;

  figure(1) % Figura S3.6
  plot(t,ysin,t,yl,[0 10],[0 0],'k');
  grid;
  xlabel('t')
  ylabel('y_L , y')
  title('Comportamento in equilibrio')
  legend('Modello non lineare','Modello linearizzato')

endfunction

function [xp] = sistema(t, x, A, w, ubar)
  % Evoluzione del sistema secondo le equazioni date
  u = ubar + A*cos(w.*t);
  xp = x.^2 - u*x - 2*u;
endfunction

