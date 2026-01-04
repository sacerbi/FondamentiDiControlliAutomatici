function Esempio2_14()
  % Si consideri il sistema dinamico con equazione di stato
  %   xdot(t) = x(t)^2 + x(t) + u(t)
  % Si studino i punti di equilibrio al variare di x
  x = -2:0.01:2;
  u1 = 1/8;
  u2 = 1/4;
  u3 = 1/2;

  y1 = x.^2 + x + u1*ones(size(x));
  y2 = x.^2 + x + u2*ones(size(x));
  y3 = x.^2 + x + u3*ones(size(x));

  plot(x,y1,x,y2,x,y3)
  hold on
  axis([-2 2 -1 2])
  plot([-2 2], [0 0], 'k')
  legend('u(t)=1/8','u(t)=1/4','u(t)=1/2')
  grid
endfunction
