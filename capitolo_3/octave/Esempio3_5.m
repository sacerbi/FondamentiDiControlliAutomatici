function Esempio3_5()
  close all
  clear all
  clc
  % Si studia il sistema massa molla, osservando la diagonalizzabilità del sistema
  pkg load control
  % Caso 1: h^2 > 4Mk
  t=0:0.1:50; % Tempo di simulazione lungo per analizzare il comportamento a regime
  u = zeros(size(t)); % Ipotizzo forza esterna agente pari a zero
  k=1;  % Progettuale, verranno studiate varie configurazioni
  h=3;  % Progettuale, verranno studiate varie configurazioni
  M=2;  % Progettuale, verranno studiate varie configurazioni
  x0=[2 1]';  % Scelta iniziale
  [A,B,C,D] = computeMat(k,h,M);
  % Autovalori
  s1 = -h/(2*M) + sqrt((h*h)/(4*M*M) - k/M);
  s2 = -h/(2*M) - sqrt((h*h)/(4*M*M) - k/M);
  % Matrice di trasformazione e sistema equivalente
  Td = 1/(s1-s2)*[s2 -1; -s1 1];
  Ad = Td * A * pinv(Td)
  xl = 1/(s2-s1) * [ (s2*x0(1)-x0(2))*e.^(s1.*t) - (s1*x0(1)-x0(2))*e.^(s2.*t);
                    (s1*s2*x0(1)-s1*x0(2))*e.^(s1.*t) - (s1*s2*x0(1)-s2*x0(2))*e.^(s2.*t)];
  % Confronto la soluzione matematica con la soluzione del risolutore
  comparazione(A,B,C,D,u,t,x0,xl(1,:),xl(2,:));

  % Caso 2: h^2 < 4Mk
  k2=1;
  h2=1;
  M2=2;
  [A2,B2,C2,D2] = computeMat(k2,h2,M2);
  % Autovalori
  s = -h2/(2*M2);
  w = sqrt(k2/M2 - (h2*h2)/(4*M2*M2));
  d = sqrt(s*s + w*w);
  g = acos(s/d); %asin(w/d);
  % Evoluzione libera
  xl2 = e.^(s.*t).*[-d/w .*sin(w.*t - g) .*x0(1) + 1/w .*sin(w.*t) .*x0(2);
                    -(d*d)/w .*sin(w*t) .*x0(1) + d/w .*sin(w.*t +g) .*x0(2)];
  % Confronto la soluzione matematica con la soluzione del risolutore
  comparazione(A2,B2,C2,D2,u,t,x0,xl2(1,:),xl2(2,:));

  % Caso 3: h^2 == 4Mk
  k3 = 2;
  h3 = 4;
  M3 = 2;
  [A3,B3,C3,D3] = computeMat(k3,h3,M3);
  % Autovalori
  s0 = -h3/(2*M3)
  % Matrice di trasformazione e sistema equivalente
  Tj = [0 -1/(s0*s0); 1 -1/s0];
  Aj = Tj * A3 * pinv(Tj)
  % Evoluzione libera
  xl3 = e.^(s0.*t) .* [x0(1) - (s0*x0(1) - x0(2)).*t; x0(2) - s0*(s0*x0(1) -x0(2)).*t];
  % Confronto la soluzione matematica con la soluzione del risolutore
  comparazione(A3,B3,C3,D3,u,t,x0,xl3(1,:),xl3(2,:));
endfunction

function [A,B,C,D] = computeMat(k,h,M)
  A=[0 1;-k/M -h/M];
  B=[0;1/M];
  C=[1 0];
  D=0;
endfunction

function comparazione(A,B,C,D,u,t,x0,x1,x2)
  % Simulo per confrontare i risultati
  sys = ss(A,B,C,D);
  [y,T,X] = lsim(sys,u,t,x0);

  %
  % plot
  %
  figure
  subplot(211),plot(t,X(:,1),';x_l(t) simulato;',t,x1,';x_l(t) teorico;');
  title('movimento libero - posizione')
  hold on
  subplot(212),plot(t,X(:,2),';dx_l(t)/dt simulato;',t,x2,';dx_l(t)/dt teorico;')
  title('movimento libero - velocità')
  xlabel('tempo')
endfunction


