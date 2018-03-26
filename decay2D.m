function sim=decay2D(v,x);
%
%
sim=v(2)*exp(-v(1)*x.^(2/3));
