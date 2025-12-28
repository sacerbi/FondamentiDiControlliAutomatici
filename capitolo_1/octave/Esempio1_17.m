function Esempio1_17()
  close all
  clear all
  clc
  % Si voglia implementare uno script che secondo l'esempio riporduca i grafici
  % della Figura 1.15a e 1.15b

  %
  % Dati 1.15: Processo massa-molla nominale sotto l'azione del controllore:
  %   Fm(t) = Fmc''(t) = mu * e(t) + nu * int_0^t (e(tau)dtau) , mu>0, nu>0
  % Si rappresentino:
  %   a) transitorio della posizione
  %   b) transitorio della forza motrice
  %
  Mbar = 1; %Volore scelto di progetto
  hbar = 3; %Volore scelto di progetto
  kbar = 1; %Volore scelto di progetto
  sbar = 1; %Volore scelto di progetto
  Fe = 0; %Ipotesi in cui Fe = Febar
  s0 = 0; %Condizioni iniziali nulle sulla posizione
  s0dot = 0; %Condizioni iniziali nulle sulla velocità
  int0 = 0; %Condizioni iniziali nulle sull'integrale d'errore
  x0 = [s0 s0dot int0]; %Vettore condizioni iniziali
  dt = 0.001; %Definizione del tempo
  t = 0:dt:10;
  mus = [1 1.5 3]; % Valori di mu per fare il grafico
  nus = [0.5 1]; % Valori di nu per fare il grafico

  % Figura 1.15a
  transitorioPosizione = zeros([length(t), length(mus)*length(nus)]);
  transitorioFm = zeros([length(t), length(mus)*length(nus)]);
  idx = 1;
  legends = {};
  for i = 1 : length(mus)
    for j = 1 : length(nus)
      [t,x] = ode45(@controlloreAnelloChiuso, t, [x0], Mbar, hbar, kbar, Fe, sbar, mus(i), nus(j));
      transitorioPosizione(:,idx) = x(:,1);
      transitorioFm (:,idx) = mus(i) * sbar - mus(i) * x(:,1) + nus(j) * x(:,3);
      legends{idx} = ['\mu = ' num2str(mus(i)) ', \nu = ' num2str(nus(j))];
      idx = idx + 1;
    endfor
  endfor

  figure(1)
  plot(t,transitorioPosizione, [0 10], [1 1], 'k');
  grid;
  xlabel('Tempo (t)')
  ylabel('Stato s')
  title('Evoluzione dello stato s per diversi guadagni (\mu, \nu) e C.I. nulle')
  legend(legends,'Location','SouthEast');

  figure(2)
  plot(t,transitorioFm, [0 10], [0 0], 'k');
  grid;
  xlabel('Tempo (t)')
  ylabel('Forza motrice')
  title('Evoluzione della Fm per diversi guadagni (\mu, \nu) e C.I. nulle')
  legend(legends,'Location','SouthEast');

endfunction

function [xp] = controlloreAnelloChiuso(t, x, M, h, k, Fe, s, mu, nu)
  % Questa funzione implementa l'evoluzione di stato per l'esempio descritto.
  % Data l'equazione
  %   M * sddot(t) = -h * sdot(t) - k * s(t) + mu * e(t) +
  %                   nu * int_0^t [e(tau)dtau] - Fe
  % Si considera inoltre e = sbar - s e lo stato composto da
  %   x(t) = [ s(t)  sdot(t) int_0^t]
  % Pertanto l'evoluzione di stato sarà
  %   xdot(t) = [ sdot(t)
  %               ((-h/M)*sdot(t) - (k/M)*s(t) + (mu/M)*e(t) + (nu/M)*int_0^t - Fe/M
  %               e(t)]
  e = s - x(1);
  xp = [x(2) ...
        -(h/M) * x(2) - (k/M) * x(1) + (mu/M) * e + (nu/M) * x(3) - Fe/M ...
        e];
endfunction
