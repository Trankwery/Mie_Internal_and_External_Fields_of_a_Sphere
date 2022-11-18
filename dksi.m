function res = dksi(l,x)
 res = -0.5*sqrt(2*pi/x).*(besselj(l+0.5,x).*l-besselj(l-0.5,x).*x+1i*bessely(l+0.5,x).*l-1i*bessely(l-0.5,x)*x);      