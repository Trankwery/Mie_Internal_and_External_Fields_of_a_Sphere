function res = fpsi(l,x)
  res = x.*( sqrt(pi/x/2)*besselj( (l+0.5),x) );
