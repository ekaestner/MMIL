function Y = rsi_SH_matrix(U,order)
%function Y = rsi_SH_matrix(U,order)
%
% Purpose: generate spherical harmonic parameterization matrix
%
% Required Input:
%   U: [nr x 3] matrix of desired reconstruction points on unit sphere 
%   order: spherical harmonic order
%
% Output:
%   Y: SH parameterization matrix [nr x nb]
%
% Created:  05/08/12 by Don Hagler
% Last Mod: 10/27/12 by Don Hagler
%

% based on code by Nate White created 09/19/07

if ~mmil_check_nargs(nargin,1), return; end;

if ~exist('order','var') || isempty(order), order = 8; end
if size(U,2)>size(U,1), U = U';end

nr = size(U,1);
nb = (order^2 + order + 2)/2 + order;

[Phi,Theta,Rho] = cart2sph(U(:,1),U(:,2),U(:,3));
Phi(find(Phi<0)) = Phi(find(Phi<0))+2*pi;
Theta = pi/2-Theta;

Y = zeros(nr,nb);
for l = 0:2:order
  p = legendre(l,cos(Theta)');
  for m = -l:l
    j = (l^2 + l + 2)/2 + m;
    if m<0
      c1 = factorial(l-m)/factorial(l+m);
      c2 = (2*l+1)/(4*pi);
      c3 = (-1)^(abs(m))*(factorial(l-abs(m))/factorial(l+abs(m)));
      Ylm = sqrt(c1*c2)*c3*p(abs(m)+1,:)'.*exp(i*m*Phi);
      Y(:,j) = sqrt(2)*real(Ylm);
    elseif m>0
      c1 = factorial(l-m)/factorial(l+m);
      c2 = (2*l+1)/(4*pi);
      Ylm = sqrt(c1*c2)*p(m+1,:)'.*exp(i*m*Phi);
      Y(:,j) = sqrt(2)*imag(Ylm);
    else
      c1 = factorial(l-m)/factorial(l+m);
      c2 = (2*l+1)/(4*pi);
      Ylm = sqrt(c1*c2)*p(m+1,:)'.*exp(i*m*Phi);
      Y(:,j) = Ylm;
    end
  end
end

return;

