function Esempio3_25()
  % Esercizio sulla risposta del sistema
  close all
  clear all
  clc
  pkg load control
  % Dati
  disp('sistema:')
  A=[-1 0;0 -2]
  B=[5;4]
  C=[1 3]
  D=0;

  sys = ss(A,B,C,D);

  % Verifica di stabilità asintotica:
  stabile = true;
  p = pole(sys);
  for i = 1:length(p)
    stabile = stabile && (p(i) < 0);
  endfor

  if(stabile) disp("sistema stabile");
  else disp("sistema NON stabile"); endif

  % Calcolo della risposta all'impulso
  figure();
  [yimp, t, x] = impulse(sys);
  plot(yimp, t);
  title ("Risposta all'impulso");

  % Calcolo della risposta alla rampa
  figure();
  [yramp, t, x] = ramp(sys);
  plot(yramp, t);
  title ("Risposta alla rampa");

  % Calcolo del guadagno statico
  % Basta guardare per s->0
  fdt = tf(sys);
  [num, den] = tfdata(fdt, 'v');
  gain = polyval(num, 0) / polyval(den, 0);
  disp("Guadagno statico:");
  gain
endfunction

