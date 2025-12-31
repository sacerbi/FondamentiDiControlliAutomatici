function Esercizio1_6()
  close all
  clear all
  clc
  % Si applichi il controllore a relè con isteresi per il quale e = hbar - h(t).
  % Si faccia riferimento al serbatoio dell'esempio 1.13.
  % Si studi l'evoluzione h(t) per h(0)=h0 generico al variare di Ehat, qhat.

  %
  % Dati 1.17: Processo di riempimento ad isteresi con portata costante qhat:
  %   A * hdot(t) = qhat , e>Ehat
  %               = -qhat , e<Ehat
  % Nel transitorio iniziale si fa riferimento all'andamento della variazione
  % di livello
  %
  A = 0.5; % Area iniziale
  qbar = 2; % Portata discreta
  Ebar = 0.5; % Margine d'errore per isteresi
  hbar = 3; % Target di livello
  h0 = 0;   % Livello iniziale
  qprev = qbar; %Variabile di appoggio per gestione isteresi portata\

  dt = 0.001; %Definizione del tempo
  ti = 0;     %Tempo inizio
  tf = 3;     %Tempo fine
  t = ti:dt:tf;  %Vettore tempi
  q = zeros(size(t)); %Vettore portata
  h = ones(size(t)) * h0; %Vettore livello

  for i = 1 : length(t)-1
    e = hbar - h(i);
    if (e >= Ebar)
      q(i) = qbar;
    elseif (e <= -Ebar)
      q(i) = -qbar;
    else
      q(i) = qprev;
    endif
    qprev = q(i);
    hdot = q(i)/A;
    h(i+1) = h(i) + hdot * dt;
  endfor

  figure(1)
  plot(t,h, [ti tf], [hbar hbar], 'k');
  ylim([0 4]);
  grid;
  xlabel('Tempo (t)')
  ylabel('Livello h(t)')
  title('Evoluzione del livello h(t) per C.I. nulle')

  figure(2)
  plot(t,q, [ti tf], [0 0], 'k');
  ylim([-3 3]);
  grid;
  xlabel('Tempo (t)')
  ylabel('Portata q(t)')
  title('Evoluzione della portata q(t) per C.I. nulle')

endfunction
