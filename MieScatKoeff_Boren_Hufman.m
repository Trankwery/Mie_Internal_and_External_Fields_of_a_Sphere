function [a,b] = MieScatKoeff_Boren_Hufman(q,n)
%  function [eBl mBl] = MieScatKoeff(q,n)
% Funkcia oblicza wspulczynniki eBl mBl
% q = 2*pi*a / lambda

% lambda  - dlugosc fali 
% n - zespolony wspulczynnik zalamania sfery wzgledem osrodka
% n = 1.33;

fjn = @(x,l) ( sqrt(pi/x/2)*besselj( (l+0.5),x) );
fpsi = @(x,l) x.*( sqrt(pi/x/2)*besselj( (l+0.5),x) );
fhi  = @(x,l) ( sqrt(pi/x/2)*bessely( (l+0.5),x) );
fksi = @(x,l) x.*( ( sqrt(pi/x/2)*besselj( (l+0.5),x) ) + complex(0,1)*( sqrt(pi/x/2)*bessely( (l+0.5),x) ) );



dpsi = @(x,l) -0.5*sqrt(2*pi/x).*(besselj(l+0.5,x).*l-besselj(l-0.5,x)*x);
dhi = @(x,l) -0.5*sqrt(2*pi/x).*(bessely(l+0.5,x)+ bessely(l+0.5,x).*l-bessely(l-0.5,x).*x);
dksi = @(x,l) -0.5*sqrt(2*pi/x).*(besselj(l+0.5,x).*l-besselj(l-0.5,x).*x+1i*bessely(l+0.5,x).*l-1i*bessely(l-0.5,x)*x);            
              
% ==========================================================================================


l = 1:ceil ( real( q ) + 4 * real( q ).^0.33333333  );

% koef = complex(0,1).^(l+1).*( ( 2*l+1 ) ./ ( l.*(l+1) ) );
cz = n*dpsi(q,l).*fpsi(n*q,l) - fpsi(q,l).*dpsi(n*q,l);
zn = n*dksi(q,l).*fpsi(n*q,l) - fksi(q,l).*dpsi(n*q,l);
a = cz./zn;

cz = fpsi(n*q,l).*dpsi(q,l) - n*dpsi(n*q,l).*fpsi(q,l);
zn = dksi(q,l).*fpsi(n*q,l) - n*fksi(q,l).*dpsi(n*q,l);
b = cz./zn;
