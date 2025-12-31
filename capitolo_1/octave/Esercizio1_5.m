function Esercizio1_5()

  close all
  clear all
  clc
  % Si applichi al processo (1.8) il controllore descritto da
  %  Fm(t) = Fmc'''(t) = kbar * sbar + Febar + mu * e(t) + nu * edot(t)
  % E si analizzi il comportamento statico e dinamico del sistema complessivo,
  % seguendo le linee sviluppate negli Esempi 1.14, 1.15, 1.17.
  Mbar = 1; %Volore scelto di progetto
  hbar = 3; %Volore scelto di progetto
  kbar = 1; %Volore scelto di progetto
  sbar = 1; %Volore scelto di progetto
  Fe = 0;   %Ipotesi in cui Fe = Febar
  Febar = 0;%Ipotesi in cui Fe = Febar
  dFe = Febar - Fe; %Ipotesi in cui Fe = Febar
  s0 = 0; %Condizioni iniziali nulle sulla posizione
  s0dot = 0; %Condizioni iniziali nulle sulla velocità
  x0 = [s0 s0dot]; %Vettore condizioni iniziali
  dt = 0.001; %Definizione del tempo
  t = 0:dt:10;
  mus = [5 10 20]; % Valori di mu per fare il grafico
  etas = [0 1 2 3 4]; % Valori di eta per fare il grafico

  % Figura 1.15a
  transitorioPosizione = zeros([length(t), length(mus)*length(etas)]);
  transitorioFm = zeros([length(t), length(mus)*length(etas)]);
  idx = 1;
  legends = {};
  for i = 1 : length(mus)
    for j = 1 : length(etas)
      [t,x] = ode45(@controlloreAnelloChiuso, t, [x0], Mbar, hbar, kbar, dFe, sbar, mus(i), etas(j));
      transitorioPosizione(:,idx) = x(:,1);
      e = sbar - x(:,1);
      edot = -x(:,2);
      transitorioFm (:,idx) = kbar*sbar + Febar + mus(i) * e + etas(j) * edot;
      legends{idx} = ['\mu = ' num2str(mus(i)) ', \eta = ' num2str(etas(j))];
      idx = idx + 1;
    endfor
  endfor

  figure(1)
  plot(t,transitorioPosizione, [0 10], [1 1], 'k');
  grid;
  xlabel('Tempo (t)')
  ylabel('Stato s')
  title('Evoluzione dello stato s per diversi guadagni (\mu, \eta) e C.I. nulle')
  legend(legends,'Location','SouthEast');

  figure(2)
  plot(t,transitorioFm, [0 10], [0 0], 'k');
  grid;
  xlabel('Tempo (t)')
  ylabel('Forza motrice')
  title('Evoluzione della Fm per diversi guadagni (\mu, \eta) e C.I. nulle')
  legend(legends,'Location','SouthEast');

endfunction

function [xp] = controlloreAnelloChiuso(t, x, M, h, k, dFe, s, mu, eta)
  % Questa funzione implementa l'evoluzione di stato per l'esempio descritto.
  % Data l'equazione
  %   M * sddot(t) = -h * sdot(t) - k * s(t) + kbar * sbar + Febar + mu * e(t) +
  %                   eta * edot(t) - Fe
  % Si considera inoltre e = sbar - s e lo stato composto da
  %   x(t) = [ s(t)  sdot(t)]
  % Pertanto l'evoluzione di stato sarà
  %   xdot(t) = [ sdot(t)
  %               ((-(h+eta)/M)*sdot(t) - ((k+mu)/M)*s(t) + ((kbar+mu)/M + dFe]
  xp = [x(2) ...
        -((h+eta)/M) * x(2) - ((k+mu)/M) * x(1) + ((k+mu)/M) * s + dFe/M];
endfunction

