function Esempio3_3()
  close all
  clear all
  clc
  % Viene implementata la rappresentazione simbolica del sistema e la sua
  % rappresentazione equivalente
  pkg load symbolic
  syms k M1 M2
  uno = sym(1);
  zero = sym(0);

  A = [ zero  zero  uno   zero;
        zero  zero  zero  uno;
        -k/M1 k/M1  zero  zero;
        k/M2  -k/M2 zero  zero]
  B = [ zero    zero;
        zero    zero;
        uno/M1  zero;
        zero    uno/M2]
  C = [ uno zero  zero  zero]

  T = [ uno zero  zero  zero;
       -uno uno   zero  zero;
       zero zero  uno   zero;
       zero zero  -uno  uno]

  At = T * A * pinv(T)
  Bt = T * B
  Ct = C * pinv(T)


endfunction
