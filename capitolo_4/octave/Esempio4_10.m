function Esempio4_10_octave()
  close all
  clear all
  clc

  % --- Condizioni iniziali che restano nel bacino di attrazione dell'origine ---
  % Vengono simulate a lungo (T = 20)
  x10 = [-2 2 -2    2];
  x20 = [ 0 0 0.1 -0.1];

  % --- Condizioni iniziali che escono dal bacino di attrazione ---
  % La traiettoria diverge (finite escape time attorno a t ~ 2.4s),
  % quindi si simula solo fino a T = 2.2, come nello script MATLAB originale
  x10a = [2 -2];
  x20a = [1 -1];

  x1 = {};
  x2 = {};
  legends = {};

  % Simulazioni con durata T = 20
  for i = 1:length(x10)
    [~, X] = ode45(@sistema, [0 20], [x10(i) x20(i)]);
    x1{i} = X(:,1);
    x2{i} = X(:,2);
    legends{i} = ['x_1(0) = ' num2str(x10(i)) ', x_2(0) = ' num2str(x20(i))];
  endfor

  % Simulazioni con durata T = 2.2 (traiettorie divergenti)
  for i = 1:length(x10a)
    [~, X] = ode45(@sistema, [0 2.2], [x10a(i) x20a(i)]);
    x1{i+length(x10)} = X(:,1);
    x2{i+length(x10)} = X(:,2);
    legends{i+length(x10)} = ['x_1(0) = ' num2str(x10a(i)) ', x_2(0) = ' num2str(x20a(i))];
  endfor

  % --- Grafico delle traiettorie nel piano di fase ---
  figure
  grid on
  hold on
  for n = 1:length(x1)
    plot(x1{n}, x2{n})
  endfor
  axis([-4 4 -4 4])
  xlabel('x_1'), ylabel('x_2')

  % --- Linea di livello V(x) = 2.2 ---
  xx = -4:0.1:4;
  yy = -4:0.1:4;
  P = [3 -1; -1 2];
  z = zeros(length(xx), length(yy));
  for i = 1:length(xx)
    for j = 1:length(yy)
      z(i,j) = 0.5 * [xx(i) yy(j)] * P * [xx(i) yy(j)]';
    endfor
  endfor
  contour(xx, yy, z', [2.2 2.2], 'r');

  % --- Linea di livello V_punto(x) = 0 ---
  zp = zeros(length(xx), length(yy));
  for i = 1:length(xx)
    for j = 1:length(yy)
      zp(i,j) = -xx(i)^2 - yy(j)^2 - xx(i)^3*yy(j) + 2*xx(i)^2*yy(j)^2;
    endfor
  endfor
  contour(xx, yy, zp', [0 0], 'g');

  legend(legends, 'Location', 'SouthEast');
  hold off
endfunction

function [xp] = sistema(t, x)
  %   x1dot(t) = -x2(t)
  %   x2dot(t) = x1(t) - x2(t) + x1(t)^2*x2(t)
  xp = [ -x(2) ...
         x(1)*x(1)*x(2) + x(1) - x(2)];
endfunction
