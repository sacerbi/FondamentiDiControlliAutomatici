function Esercizio1_3()

  close all
  clear all
  clc

  %
  % Dati 1.3: [...] si consideri ancora l'equazione (1.10) e si studi il caso in
  % cui sia mu <= -kbar, con |mu| elevato, assumendo dapprima che il processo
  % sia in condizioni nominali e successivamente che, rimanendo fissati ai valori
  % nominali gli altri parametri, si abbia Fe != Febar
  %
  Mbar = 1; %Volore scelto di progetto
  hbar = 3; %Volore scelto di progetto
  kbar = 1; %Volore scelto di progetto
  sbar = 1; %Volore scelto di progetto
  dFe = 0; %Ipotesi in cui Fe = Febar
  s0 = 0.2; %Condizioni iniziali non nulle
  s0dot = 0.2; %Condizioni iniziali non nulle
  x0 = [s0 s0dot]; %Vettore condizioni iniziali
  dt = 0.001; %Definizione del tempo
  t = 0:dt:10; %Vettore tempo
  mus = [ -kbar, (-kbar -0.01) ];  % Primo caso. mu = -kbar in condizioni nominali
                                % Secondo caso. mu < -kbar in cond. nom.

  %% Prima richiesta

  risultatiControllore = zeros([length(t), length(mus)]);
  legends = {};
  for i = 1 : length(mus)
    [t,x] = ode45(@controlloreAnelloChiuso, t, [x0], Mbar, hbar, kbar, dFe, sbar, mus(i));
    risultatiControllore(:,i) = x(:,1);
    legends{i} = ['\mu = ' num2str(mus(i))];
  endfor

  figure(1)
  plot(t,risultatiControllore, [0 10], [1 1], 'k');
  grid;
  xlabel('Tempo (t)')
  ylabel('Stato s')
  title('Evoluzione dello stato s per diversi guadagni (\mu) e C.I. nominali')
  legend(legends,'Location','SouthEast');

  %% Seconda richiesta

  dFe = 0.2; %Ipotesi in cui Fe != Febar
  x0 = [s0 s0dot]; %Vettore condizioni iniziali
  t = 0:dt:100; %Vettore tempo

  risultatiControllore2 = zeros([length(t), length(mus)]);
  legends = {};
  for i = 1 : length(mus)
    [t,x] = ode45(@controlloreAnelloChiuso, t, [x0], Mbar, hbar, kbar, dFe, sbar, mus(i));
    risultatiControllore2(:,i) = x(:,1);
    legends{i} = ['\mu = ' num2str(mus(i))];
  endfor

  figure(2)
  plot(t,risultatiControllore2, [0 10], [0 0], 'k');
  grid;
  xlabel('Tempo (t)')
  ylabel('Stato s')
  title('Evoluzione dello stato s per diversi guadagni (\mu) e perturbazione Fe')
  legend(legends,'Location','SouthEast');

endfunction

function [xp] = controlloreAnelloChiuso(t, x, M, h, k, dFe, s, mu)
  % Questa funzione implementa l'evoluzione di stato per l'esempio descritto.
  % Data l'equazione
  %   M * sddot(t) = -h * sdot(t) - (k+mu)*s(t) + (kbar + mu)*sbar - Fe + Febar
  % Si considera Febar - Fe = dFe (delta Fe). Inoltre si assume lo stato come
  %   x(t) = [ s(t)  sdot(t) ]
  % Pertanto l'evoluzione di stato sarà
  %   xdot(t) = [ sdot(t) ((-h/M)*sdot(t) - ((k+mu)/M)*s(t) + ((kbar + mu)/M)*sbar + dFe/M)]
  xp = [x(2) ...
        -(h/M) * x(2) - ((k + mu)/M * x(1)) + ((k + mu)/M) * s + dFe/M];
endfunction

