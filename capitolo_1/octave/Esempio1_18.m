function Esempio1_18()
  close all
  clear all
  clc
  % Si voglia implementare uno script che secondo l'esempio riporduca i grafici
  % della Figura 1.17a e 1.17b

  %
  % Dati 1.17: Processo massa-molla nominale sotto l'azione del controllore ad isteresi:
  %   Fm(t) = Fm2 , e>E2
  %         = Fm1 , e<E1
  % Si rappresentino:
  %   a) transitorio della posizione
  %   b) transitorio della forza motrice
  %
  Mbar = 1; %Valore scelto di progetto
  hbar = 3; %Valore scelto di progetto
  kbar = 1; %Valore scelto di progetto
  sbar = 1; %Valore scelto di progetto
  Febar = 0; %Valore scelto di progetto
  E1=-0.02; %Valore scelto di progetto
  E2=0.04; %Valore scelto di progetto
  Fm1=-1; %Valore scelto di progetto
  Fm2=2; %Valore scelto di progetto
  Fm=0; %Valore scelto di progetto
  s0 = 0; %Condizioni iniziali nulle sulla posizione
  s0dot = 0; %Condizioni iniziali nulle sulla velocità
  x0 = [s0 s0dot]; %Vettore condizioni iniziali

  dt = 0.001; %Definizione del tempo
  t = 0:dt:10;

  % Figura 1.15a
  [t,x] = ode45(@controlloreAnelloChiuso, t, [x0], Mbar, hbar, kbar, Febar, sbar, E1, Fm1, E2, Fm2, Fm);
  transitorioPosizione = x(:,1);
  transitorioVelocita = x(:,2);
  transitorioFm = zeros(size(transitorioPosizione));

  for i = 1 : length(transitorioFm)
    errore = sbar - transitorioPosizione(i);
    if errore>=E2
      transitorioFm(i)=Fm2;
    elseif errore<=E1
      transitorioFm(i)=Fm1;
    elseif transitorioVelocita(i)>=0
      transitorioFm(i)=Fm2;
    else
      transitorioFm(i)=Fm1;
    end
  endfor

  figure(1)
  plot(t,transitorioPosizione, [0 10], [1 1], 'k');
  grid;
  xlabel('Tempo (t)')
  ylabel('Stato s')
  title('Evoluzione dello stato s per C.I. nulle')

  figure(2)
  plot(t,transitorioFm, [0 10], [0 0], 'k');
  grid;
  xlabel('Tempo (t)')
  ylabel('Forza motrice')
  title('Evoluzione della Fm per C.I. nulle')

endfunction

function [xp] = controlloreAnelloChiuso(t, x, M, h, k, Fe, s, E1, Fm1, E2, Fm2, Fm)
  % Questa funzione implementa l'evoluzione di stato per l'esempio descritto.
  % Data l'equazione
  %   M * sddot(t) = -h * sdot(t) - k * s(t) + Fm(t) - Fe
  % Si considera inoltre e = sbar - s e lo stato composto da
  %   x(t) = [ s(t)  sdot(t)]
  % Pertanto l'evoluzione di stato sarà
  %   xdot(t) = [ sdot(t)
  %               ((-h/M)*sdot(t) - (k/M)*s(t) + (1/M)*Fm(t) - Fe/M ]
  % Il valore di Fm(t) è Fm2 quando l'errore supera E2 e ci rimane finchè non
  % scende sotto E1. All'inizio viene inizializzato a seconda della velocità
  e = s - x(1);
  if e>=E2
    Fm=Fm2;
  elseif e<=E1
    Fm=Fm1;
  elseif x(2)>=0
    Fm=Fm2;
  else
    Fm=Fm1;
  end

  xp = [x(2) ...
        -(h/M) * x(2) - (k/M) * x(1) + (1/M) * Fm - Fe/M ];
endfunction
