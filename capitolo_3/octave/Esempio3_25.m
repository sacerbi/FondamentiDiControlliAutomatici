function Esempio3_25()
  % Esercizio sulla raggiungibilità
  close all
  clear all
  clc
  % Dati
  disp('sistema:')
  R=1;
  C=2;
  A=[-1/(R*C) -1/(R*C);-1/(R*C) -1/(R*C)]
  B=[1/(R*C);1/(R*C)]
  C=[1 0]
  D=0;

  % Matrice di raggiungibilità
  disp('matrice di raggiungibilità:')
  Mr = [B A*B]

  % Verifica di raggiungibilità
  nr = rank(Mr);
  if nr == size(B,1)
    disp('Completamente raggiungibile');
    return;
  else
    disp('Non completamente raggiungibile');
  endif

  % Se non completamente raggiungibile, scompongo
  disp('matrice di trasformazione:')
  T1=[Mr(:,1) [0;1]]; % Abbiamo n = 2 e rk = 1
  T=inv(T1)

  % Scompongo
  Ar = T*A*T1
  Br = T*B
  Cr = C*T1
endfunction

