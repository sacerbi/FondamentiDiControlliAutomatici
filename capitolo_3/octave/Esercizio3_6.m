function Esercizio3_6()
  close all
  clear all
  clc
  % L'obiettivo è applicare il criterio di Kharitonov
  phi_min = [3, 2, 6, 1];
  phi_max = [4, 3, 6, 2];

  order = size(phi_max,2);
  Pa = []; Pb = []; Pc = []; Pd = [];
  % sfrutto degli indici per capire quando prendere phi_max e quando phi_min
  idx = 1 : order;
  % Creo pattern di 0 e 1 secondo il pattern del criterio
  idx_d = (mod(idx, 4) >= 2);                 % [- + + - - + + - ...]
  idx_c = (mod(idx,4) < 2);                   % [+ - - + + - - + ...]
  idx_b = (mod(idx,4) == 0 | mod(idx,4) == 3);% [- - + + - - + + ...]
  idx_a = (mod(idx,4) == 1 | mod(idx,4) == 2);% [+ + - - + + - - ...]

  % Funzione che compone i polinomi. Se l'indice i-esimo passato è 0, verrà preso
  % il valore minimo, mentre se è 1 verrà preso il valore massimo.
  Pa = componi(idx_a, phi_min, phi_max);
  Pb = componi(idx_b, phi_min, phi_max);
  Pc = componi(idx_c, phi_min, phi_max);
  Pd = componi(idx_d, phi_min, phi_max);

  Pa
  Pb
  Pc
  Pd

  sa = isStabile(routh(Pa, 0.001));
  sb = isStabile(routh(Pb, 0.001));
  sc = isStabile(routh(Pc, 0.001));
  sd = isStabile(routh(Pd, 0.001));

  stabile = sa && sb && sc && sd

endfunction

function stabile = isStabile(rt)
  stabile = all(rt(:,1)>0) || all(rt(:,1)<0);
endfunction

function poli = componi(indici, minVals, maxVals)
  poli = (1-indici).*minVals + indici .* maxVals;
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


