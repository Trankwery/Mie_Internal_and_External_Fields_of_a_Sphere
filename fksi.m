function res = fksi(l,x)
 res = x.*( ( sqrt(pi/x/2)*besselj( (l+0.5),x) ) + complex(0,1)*( sqrt(pi/x/2)*bessely( (l+0.5),x) ) );