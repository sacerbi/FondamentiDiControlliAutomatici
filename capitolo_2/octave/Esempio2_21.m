function Esempio2_21()
  close all
  clear all
  clc
  % Si voglia riprodurre un grafico x1-x2 rappresentante l'evoluzione del sistema
  % identificato dal circuito dell'esempio, nelle condizioni in qui questo si
  % specializza in un oscillatore di van der Pol.
  C = 1;
  L = 1;
  alpha = 1;
  beta = 1/3;
  x1_0=[-3 -2.5 -2 3  3 2.5 2 -0.01 0.01];
  x2_0=[ 4   4   4 4 -4 -4 -4   0    0  ];
  t = 0:0.001:40;

  evoluzioniX1 = {};
  evoluzioniX2 = {};
  legends = {};
  for i = 1 : length(x1_0)
    [T1, X1] = ode45(@oscillatore, t, [x1_0(i) x2_0(i)], C, L, alpha, beta);
    evoluzioniX1{i} = X1(:,1);
    evoluzioniX2{i} = X1(:,2);
    legends{i} = ['x_1(0) = ' num2str(x1_0(i)) ', x_2(0) = ' num2str(x2_0(i))];
  endfor

  figure
  grid
  hold on
  for i = 1 : length(x1_0)
    x1 = evoluzioniX1{i};
    x2 = evoluzioniX2{i};
    plot(x1, x2)
  endfor
  axis([-4 4 -4 4])
  xlabel('x_1'),ylabel('x_2')
  legend(legends,'Location','SouthEast');

endfunction

function [xp] = oscillatore(t, x, C, L, alpha, beta)
  % Per questo oscillatore, si considera x1(t) = vc(t) e x2(t) = ic(t).
  % Dalle relazioni i = C * (dvc /dt) e vl = L * (dil/dt) è possibile,
  % in funzione della legge ic + id + il = 0 al nodo, ricavare le eq. per le
  % derivate degli stati x1 e x2:
  %   x1dot(t) = (1/C) * x2(t)
  %   x2dot(t) = -[3 * beta * x1(t)^2 - alpha] * (x2(t)/C) - 1/L * x1(t)
  xp = [ x(2)/C ...
         -(3 * beta * x(1)^2 - alpha) * (x(2)/C) - (x(1)/L)];
endfunction
