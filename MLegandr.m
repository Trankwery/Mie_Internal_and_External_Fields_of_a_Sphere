function [dP P] = MLegandr(theta,l)
%
% theta [ radian ]
%

tP = zeros(l+1,length(theta));
tdP = tP;

tP(2,:) = 1;
tdP(2,:) = cos(theta);
j = 2;

for ii = 3:(l+1)
    tP(ii,:) =( ( cos(theta)*(2*j-1)/(j-1) ).* tP(ii-1,:) ) - ( ( j/(j - 1) )* tP(ii - 2, :) );
    tdP(ii,:) = ( cos(theta)*j .* tP( ii ,:) ) -  ( (j + 1) .* tP(ii-1,:));
    j = j + 1;
end;

P = tP(2:end,:);
dP = tdP(2:end,:);

