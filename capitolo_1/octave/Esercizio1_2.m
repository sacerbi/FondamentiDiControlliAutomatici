% Esercizio 1.2 %
function Esercizio1_2()

  close all
  clear all
  clc
  % Valutando la condizione statica, la legge Cdot(t) = -aC(t) + bI(t) diventa
  %                     0 = -aCbar + bIbar
  % Il sistema ad anello aperto è quindi descritto dalla relazione
  %           Cdot(t) = -aC(t) + (abar/bbar)*bbar*Cbar
  % Risulta chiaro che C = Cbar solamente se a e b assumono i valori nominali abar
  % e bbar. Inoltre le variazioni di a influenzano i transitori:
  % se si assume la condizione iniziale C(0) = C0, la soluzione per integrazione:
  % C(t) = e^(-at)C0 + (abar*b)/(a*bbar)*Cbar(1 - e^(-at))
  %
  % Per valutare gli effetti ad anello chiuso dei controllori proporzionale e
  % integrale, bisogna risolvere il sistema usando il metodo ode45.

  % Dati %
  a=1;
  b=1;
  C0=0.2;
  I0=0.1;
  Cbar=5;
  %
  % Controllore proporzionale %
  %
  t=0:0.01:2;
  x0=0;
  mus = [1, 2, 4, 8, 16, 32];
  Cprop = zeros([length(t), length(mus)]);
  legends = {};
  for i = 1 : length(mus)
    [t,x] = ode45(@controlloreProporzionale, t, [x0], a, b, mus(i), Cbar);
    Cprop(:,i) = x(:);
    legends{i} = ['\mu = ' num2str(mus(i))];
  endfor

  figure(1)
  plot(t,Cprop);
  grid;
  xlabel('Tempo (t)')
  ylabel('Stato C')
  title('Evoluzione dello stato C per diversi guadagni (\mu)')
  legend(legends,'Location','SouthEast');

  %
  % Controllore integrale %
  %
  t=0:0.01:12;
  x0=[0 0];
  mus = mus/10;
  Cint = zeros([length(t), length(mus)]);
  legendsIntegreal = {};
  for i = 1 : length(mus)
    [t,x] = ode45(@controlloreIntegrale, t, [x0], a, b, mus(i), Cbar);
    Cint(:,i) = x(:,1);
    legendsIntegreal{i} = ['\mu = ' num2str(mus(i))];
  endfor

  figure(2)
  plot(t,Cint,[0 12], [5 5], 'k');
  grid;
  xlabel('Tempo (t)')
  ylabel('Stato C')
  title('Evoluzione dello stato C per diversi guadagni (\mu)')
  legend(legendsIntegreal,'Location','SouthEast');
endfunction

function [xp] = controlloreProporzionale(t, x, a, b, mu, Cbar)
  % In questo caso lo stato è solo C(t). Sostituendo I(t)=mu(Cbar - C(t))
  % nell'equazione di Cdot(t) = -aC(t) + bI(t) si ottiene:
  %     Cdot(t) = -(a + b*mu) * C(t) + (b * mu * Cbar)
  xp = -(a + b*mu) * x + (b * mu * Cbar);
endfunction

function [xp] = controlloreIntegrale(t, x, a, b, mu, Cbar)
  % In questo caso lo stato è composto da C(t) e I(t). Derivando la legge di consumo
  % e inserendo la legge del controllore Idot(t) = mu*(Cbar - C(t)) si ottiene:
  %     Cddot(t) + a*Cdot(t) + b*mu*C(t) = b*mu*Cbar
  % Lo stato sarà x = [C(t); I(t)]
  xp = [ -a * x(1) + b * x(2); % Coefficienti per evoluzione di C(t)
         mu * Cbar - mu * x(1)]; % Coefficienti per evoluzione di I(t)
endfunction
