function Esempio4_6()
  %
  % Disegna le traiettorie e le linee di livello della funzione di Lyapunov
  % per l'esempio 4.6
  %
  clear all
  close all
  clc
  % Sistema di equazioni differenziali per x_1 e x_2
  ode_system = @(t, x) [
      -x(1)^3 + x(1)*x(2)^2 - x(2)^2;
       x(1)*x(2) - x(1)^2*x(2) - x(2)^3
  ];

  % Vettori degli stati iniziali (presi dal grafico riportato a pag. 103)
  x10 = [-3 -3 -3 -3 -3 -1.5 -1.5 -1 -1  0 0  1 1  2 2  3  3 3 3 3];
  x20 = [-3 -1  0  1  3 -3    3   -3  3 -3 3 -3 3 -3 3 -3 -1 0 1 3];

  tspan = [0 25]; % Intervallo di simulazione

  % Simulazione e calcolo delle traiettorie
  for i = 1:length(x10)
      x0 = [x10(i); x20(i)];  % Condizione iniziale corrente
      [t, x_out] = ode45(ode_system, tspan, x0); % Simulazione
      % Salvataggio dei risultati in array di celle
      x1{i} = x_out(:,1)';
      x2{i} = x_out(:,2)';
  end

  % Visualizzazione di tutte le traiettorie nello stesso grafico
  figure;
  grid on;
  hold on;

  for n = 1:length(x1)
     plot(x1{n}, x2{n}, 'b'); % Traiettorie in blu per leggibilità
     axis([-3 3 -3 3]);
     xlabel('x_1');
     ylabel('x_2');
     title('Traiettorie e linee di livello (Octave)');
  end

  % Plotto anche le linee di livello della funzione di Lyapunov V = 0.5*(x1^2 + x2^2)
  punti = -3:0.01:3;
  [X_grid, Y_grid] = meshgrid(punti, punti);
  Z = 0.5 * (X_grid.^2 + Y_grid.^2);

  % Tracciamento curve di livello in rosso
  contour(X_grid, Y_grid, Z, [0.125 0.25 0.5 1 2 4 8], 'r');
  axis square;
  hold off;

endfunction
