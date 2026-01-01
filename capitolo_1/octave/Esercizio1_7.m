function Esercizio1_7()
  % Si consideri il processo a tempo discreto con variabile controllata y e
  % variabile di controllo u:
  %   y(k+1) = u(k)
  % e si applichi il controllore
  %   u(k) = mu * (w - y(k))
  % dove w è il target e mu!=1 è oggetto di scelta. Determinare l'equazione alle
  % differenze del sistema, ricavare il valore alle condizioni statiche per il
  % valore ybar. Si determini la generalizzazione dell'evoluzione di y(k) e si
  % determini il valore di convergenza.

  % Dati
  y0 = 0;                       % Valore iniziale del sistema
  k = 4;                        % Istante discreto fino a cui calcolare
  ws = [1 2 3];                     % Target
  mus = [-1.1, -0.5, 0.5, 1.1]; % Valori del controllore
  t = 0:k;                      % Vettore dei tempi

  % Calcolo le evoluzioni e le stampo
  for i = 1:length(ws)
    w = ws(i);
    evoluzioneProcesso = zeros([length(t), length(mus)]);
    legends = {};
    for j = 1:length(mus)
      y = evoluzioneDiscreta(y0, w, mus(j), k);
      evoluzioneProcesso(:, j) = y;
      legends{j} = ['\mu = ' num2str(mus(j)) ', w = ' num2str(w)];
    endfor
    figure(i);
    plot(t,evoluzioneProcesso, '-*-', [0 k], [0 0], 'k',
          [0 k], [w w], '--', [0 k], [-w -w], '--');
    grid;
    xlabel('Passi (k)')
    ylabel('Evoluzione y')
    title('Evoluzione del processo y per diversi guadagni (\mu) e C.I. nulle')
    legend(legends,'Location','SouthEast');
  endfor
endfunction

function y = evoluzioneDiscreta(y0, w, mu, k)
  % Calcola l'evoluzione del sistema discreto con valore iniziale y0, target w,
  % guadagno del controllore mu e numero di passi discreti k.
  % Il processo è descritto da
  %   y(k+1) = mu * (w - y(k))
  if (k == 0)
    y = y0;
    return;
  endif

  y = zeros([k+1, 1]);% Vettore dell'evoluzione
  y(1) = y0;          % Impongo la condizione iniziale

  for i = 1 : k
    y(i+1) = mu * (w - y(i));
  endfor
endfunction
