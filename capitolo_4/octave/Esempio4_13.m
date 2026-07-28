function Esempio4_13()
  close all
  clear all
  clc

  % --- Condizioni iniziali che restano nel bacino di attrazione dell'origine ---
  % Vengono simulate a lungo (T = 20)
  x10 = [0.1 0.4 0.8 1.6 2.3 2.6   4   4   4   4   4   4];
  x20 = [  4   4   4   4   4   4 3.8 3.1 2.2 1.3 0.7 0.2];


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

  % --- Grafico delle traiettorie nel piano di fase ---
  figure
  grid on
  hold on
  for n = 1:length(x1)
    plot(x1{n}, x2{n})
  endfor
  axis([0 4 0 4])
  xlabel('x_1'), ylabel('x_2')

  % --- Linea di livello V(x) ---
  xx = 0:0.01:4;
  yy = 0:0.01:4;
  z = zeros(length(xx), length(yy));
  for i = 1:length(xx)
    for j = 1:length(yy)
      z(i,j) = xx(i) - 2*(1 + log(xx(i)/2)) + yy(j) - 1*(1 + log(yy(j)/1));
    endfor
  endfor
  contour(xx, yy, z', [0.1 0.5 1 2 3], 'g');

  legend(legends, 'Location', 'SouthEast');
  hold off
endfunction

function [xp] = sistema(t, x)
  %   x1dot(t) = x1(t)(3 - x1(t) - x2(t))
  %   x2dot(t) = x2(t)(x1(t) - x2(t) -1)
  xp = [ x(1)*(3 - x(1) - x(2)) ...
         x(2)*(x(1) - x(2) - 1)];
endfunction
