function Esempio3_16()
  close all
  clear all
  clc
  % L'obiettivo è rappresentare la regione di stabilità asintotica per un
  % sistema con parametri incerti.
  alpha = -2:0.1:4;
  beta = -2:0.1:4;
  statiStabili = zeros(length(alpha),length(beta));

  for i = 1:length(alpha)
    for j = 1:length(beta)
      rt = routh([1, 2+beta(j), 1+2*beta(j), alpha(i)+beta(j)], 0.001);
      statiStabili(i,j) = all(rt(:,1)>0) || all(rt(:,1)<0);
    endfor
  endfor

  [B, A] = meshgrid(beta, alpha);

  figure(1)
  plot(A(statiStabili==1),B(statiStabili==1), 'r*',...
        A(statiStabili==0),B(statiStabili==0), 'b+');
  xlabel('alpha')
  ylabel('beta')
  title('Regione di stabilità (in rosso)')

endfunction

function rt = routh(poli,epsilon)
  % Compongo la tabella di Routh
  if (nargin < 2)
    fprintf('\nEpsilon not provided. Assuming 0.0001.');
    epsilon = 0.0001;
  endif

  nCoeff = size(poli,2);
  rt = zeros(nCoeff, ceil(nCoeff/2)); % Inizializzo la tabella

  % Si popolano le prime due righe
  for i = 1 : nCoeff
    r = 2 - mod(i,2); %Prima o seconda riga
    c = ceil(i/2); %Procedo con le colonne
    rt(r, c) = poli(i);
  endfor

  righeDaFare = nCoeff - 2;
  indiceMaxColonna = zeros(righeDaFare,1);
  for i = 1 : righeDaFare
    indiceMaxColonna(righeDaFare-i+1) = ceil(i/2);
  endfor

  for i = 3 : nCoeff
    % Prima di calcolare il nuovo valore, verifico se ci sono casi particolari.
    % I casi particolari sono descritti in:
    % https://it.wikipedia.org/wiki/Criterio_di_Routh-Hurwitz
    if (all(rt(i-1,:)==0)) % Caso di riga nulla precedente
      fprintf('\nCaso speciale riga nulla\n');
      ordPoliAux = nCoeff - i + 2; % Ordine del polinomio ausiliario
      nCoeffAux = ceil(ordPoliAux/2) - mod(ordPoliAux,2) + 1; % Num. coeff. polinomio ausiliario
      poliAux = rt(i-2,1:nCoeffAux); % Prendo il polinomio ausiliario da usare
      potenzeAux = ordPoliAux:-2:0; % Il polinomio ha potenze pari
      rt(i-1,1:nCoeffAux) = poliAux .* potenzeAux; % Assegno la derivata
    elseif (rt(i-1,1)==0) % Caso di primo termine nullo
      fprintf('\nCaso speciale primo termine nullo\n');
      if (all(rt(1:i-1,1)>=0))
        rt(i-1,1) = epsilon;
      else
        rt(i-1,1) = -epsilon;
      endif
    endif

    for j = 1 : indiceMaxColonna(i-2)
      h1 = rt(i-2,1);
      k1 = rt(i-1,1);
      hi1 = rt(i-2,j+1);
      ki1 = rt(i-1,j+1);
      rt(i,j) = hi1 - (h1*ki1)/k1; % Applico la formula di p.69
    endfor
  endfor
endfunction


