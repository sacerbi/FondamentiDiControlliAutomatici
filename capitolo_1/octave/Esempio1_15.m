function Esempio1_15()

  close all
  clear all
  clc
  % Si voglia implementare uno script che secondo l'esempio riporduca i grafici
  % della Figura 1.12, 1.13 (a e b) e 1.14

  %
  % Dati 1.12: Transitori della posizione del processo massa-molla nominale con
  % condizioni iniziali nulle sotto l'azione del controllore. Sono noti i valori
  % nominali Mbar, kbar e hbar (Esempio 1.14)
  %
  Mbar = 1; %Volore scelto di progetto
  hbar = 3; %Volore scelto di progetto
  kbar = 1; %Volore scelto di progetto
  sbar = 1; %Volore scelto di progetto
  dFe = 0; %Ipotesi in cui Fe = Febar
  s0 = 0; %Condizioni iniziali nulle
  s0dot = 0; %Condizioni iniziali nulle
  x0 = [s0 s0dot]; %Vettore condizioni iniziali
  dt = 0.001; %Definizione del tempo
  mus = [0.001 0.5 1 5 10 20]; % Valori di mu per fare il grafico

  % Figura 1.12
  t = 0:dt:10; %Vettore tempo
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
  title('Evoluzione dello stato s per diversi guadagni (\mu) e C.I. nulle')
  legend(legends,'Location','SouthEast');

  %
  % Dati 1.13a: Transitori della posizione del processo massa-molla nominale con
  % condizioni iniziali non nulle per s0 e s0dot sotto l'azione del controllore.
  %
  Mbar = 1; %Volore scelto di progetto
  hbar = 3; %Volore scelto di progetto
  kbar = 1; %Volore scelto di progetto
  sbar = 0; %Volore scelto di progetto
  dFe = 0; %Ipotesi in cui Fe = Febar
  dt = 0.001; %Definizione del tempo
  s0 = 1; %Condizioni iniziali non nulle
  s0dot = 1; %Condizioni iniziali non nulle
  x0 = [s0 s0dot]; %Vettore condizioni iniziali

  % Figura 1.13a
  t = 0:dt:10; %Vettore tempo
  mus = [0.001 0.5 1 5 10 20]; % Valori di mu per fare il grafico
  risultatiControllore = zeros([length(t), length(mus)]);
  legends = {};
  for i = 1 : length(mus)
    [t,x] = ode45(@controlloreAnelloChiuso, t, [x0], Mbar, hbar, kbar, dFe, sbar, mus(i));
    risultatiControllore(:,i) = x(:,1);
    legends{i} = ['\mu = ' num2str(mus(i))];
  endfor

  figure(2)
  plot(t,risultatiControllore, [0 10], [0 0], 'k');
  grid;
  xlabel('Tempo (t)')
  ylabel('Stato s')
  title('Evoluzione dello stato s per diversi guadagni (\mu) e C.I. non nulle')
  legend(legends,'Location','SouthEast');

  %
  % Dati 1.13b: Transitori della posizione del processo massa-molla nominale con
  % condizioni iniziali nulle per s0 e s0dot, Fe!=Febar sotto l'azione del controllore.
  %
  Mbar = 1; %Volore scelto di progetto
  hbar = 3; %Volore scelto di progetto
  kbar = 1; %Volore scelto di progetto
  sbar = 0; %Volore scelto di progetto
  dFe = 1; %Ipotesi in cui Fe = Febar
  dt = 0.001; %Definizione del tempo
  s0 = 0; %Condizioni iniziali non nulle
  s0dot = 0; %Condizioni iniziali non nulle
  x0 = [s0 s0dot]; %Vettore condizioni iniziali

  % Figura 1.13b
  t = 0:dt:10; %Vettore tempo
  mus = [0.001 0.5 1 5 10 20]; % Valori di mu per fare il grafico
  risultatiControllore = zeros([length(t), length(mus)]);
  legends = {};
  for i = 1 : length(mus)
    [t,x] = ode45(@controlloreAnelloChiuso, t, [x0], Mbar, hbar, kbar, dFe, sbar, mus(i));
    risultatiControllore(:,i) = x(:,1);
    legends{i} = ['\mu = ' num2str(mus(i))];
  endfor

  figure(3)
  plot(t,risultatiControllore, [0 10], [0 0], 'k');
  grid;
  xlabel('Tempo (t)')
  ylabel('Stato s')
  title('Evoluzione dello stato s per diversi guadagni (\mu) e C.I. non nulle')
  legend(legends,'Location','SouthEast');

  %
  % Dati 1.14: Transitori della forza motrice del processo massa-molla nominale con
  % condizioni iniziali nulle per s0 e s0dot sotto l'azione del controllore.
  % Dall'esempio 1.14 risulta che
  %   Fm = kbar * sbar + Febar + mu * e
  % dove e = sbar - s
  %
  Mbar = 1; %Volore scelto di progetto
  hbar = 3; %Volore scelto di progetto
  kbar = 1; %Volore scelto di progetto
  sbar = 1; %Volore scelto di progetto
  dFe = 0; %Ipotesi in cui Fe = Febar
  dt = 0.001; %Definizione del tempo
  s0 = 0; %Condizioni iniziali non nulle
  s0dot = 0; %Condizioni iniziali non nulle
  x0 = [s0 s0dot]; %Vettore condizioni iniziali

  % Figura 1.14
  t = 0:dt:10; %Vettore tempo
  mus = [0.001 0.5 1 5 10 20]; % Valori di mu per fare il grafico
  transitoriFm = zeros([length(t), length(mus)]);
  legends = {};
  for i = 1 : length(mus)
    [t,x] = ode45(@controlloreAnelloChiuso, t, [x0], Mbar, hbar, kbar, dFe, sbar, mus(i));
    s = x(:,1);
    transitoriFm(:,i) = kbar * sbar + dFe + mus(i) * (sbar - s);
    legends{i} = ['\mu = ' num2str(mus(i))];
  endfor

  figure(4)
  plot(t,transitoriFm, [0 10], [0 0], 'k');
  grid;
  xlabel('Tempo (t)')
  ylabel('Forza motrice')
  title('Transitorio Fm per diversi guadagni (\mu) e C.I. nulle')
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

