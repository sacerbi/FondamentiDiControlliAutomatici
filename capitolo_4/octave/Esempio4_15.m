function Esempio4_10()
  close all
  clear all
  clc

  % --- Condizioni iniziali che restano nel bacino di attrazione dell'origine ---
  % Vengono simulate a lungo (T = 20)
  x10 = [-0.99, -2.0, -2.0, -1.5, -2.0, -2.0, -1.5, -1.0, -1.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0];
  x20 = [0.0, 1.48, 1.8, 2.0, -1.48, -1.8, -2.0, 2.0, -2.0, 2.0, -2.0, 2.0, -2.0, 1.5, -1.5];

  % --- Condizioni iniziali che escono dal bacino di attrazione ---
  % La traiettoria diverge (finite escape time attorno a t ~ 4s),
  % quindi si simula solo fino a T = 4
  x10a = [-2.0, -2.0, -2.0, -2.0, -1.01];
  x20a = [-1.45, 1.45, -1.42, 1.42, 0.0];

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

  % Simulazioni con durata T = 4 (traiettorie divergenti)
  for i = 1:length(x10a)
    [~, X] = ode45(@sistema, [0 4], [x10a(i) x20a(i)]);
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
  axis([-2 2 -2 2])
  xlabel('x_1'), ylabel('x_2')

  % --- Linea di livello V(x) ---
  xx = -2:0.1:2;
  yy = -2:0.1:2;
  z = zeros(length(xx), length(yy));
  for i = 1:length(xx)
    for j = 1:length(yy)
      z(i,j) = 0.5 * [xx(i)*xx(i) + yy(j)*yy(j)];
    endfor
  endfor
  contour(xx, yy, z', [0.25,0.5,2], 'r');

  legend(legends, 'Location', 'SouthEastOutside');
  hold off
endfunction

function [xp] = sistema(t, x)
  %   x1dot(t) = -x1(t)(x1(t) + u(t)) + x2(t)^2
  %   x2dot(t) = -(x1(t) + 2u(t))x2(t)
  xp = [ -x(1)*(x(1) + 1) + x(2)*x(2) ...
        -(x(1) + 2*1)*x(2)];
endfunction
