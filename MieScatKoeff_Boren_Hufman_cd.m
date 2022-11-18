function [c,d] = MieScatKoeff_Boren_Hufman_cd(q,n)
%  function [eBl mBl] = MieScatKoeff(q,n)
% Funkcia oblicza wspulczynniki eBl mBl
% q = 2*pi*a / lambda

% lambda  - dlugosc fali 
% n - zespolony wspulczynnik zalamania sfery wzgledem osrodka
% n = 1.33;

fpsi = @(x,l) x.* sqrt(pi/(x*2))*besselj( (l+0.5),x);
fksi = @(x,l) x.* sqrt(pi/(x*2))*( besselj( (l+0.5),x)  + complex(0,1)*bessely( (l+0.5),x) );

dpsi = @(x,l) -sqrt(pi/(x*2)).*(besselj(l+1/2,x).*l-besselj(l-1/2,x)*x);
dksi = @(x,l) -sqrt(pi/(x*2)).*(besselj(l+1/2,x).*l-besselj(l-1/2,x).*x+i*bessely(l+1/2,x).*l-i*bessely(l-1/2,x)*x);            

              
% ==========================================================================================


l = 1:ceil ( real( q ) + 4 * real( q ).^0.33333333  );

% koef = complex(0,1).^(l+1).*( ( 2*l+1 ) ./ ( l.*(l+1) ) );
cz = n*fpsi(q,l).*dksi(q,l) - n*fksi(q,l).*dpsi(q,l);
zn = dksi(q,l).*fpsi(n*q,l) - n*fksi(q,l).*dpsi(n*q,l);
c = cz./zn;

cz = n*fpsi(q,l).*dksi(q,l) - n*fksi(q,l).*dpsi(q,l);
zn = n*fpsi(n*q,l).*dksi(q,l) - fksi(q,l).*dpsi(n*q,l);
d = cz./zn;