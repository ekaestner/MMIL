function R = rsi_FOD_matrix(Q,X,lADC,tADC)
%function R = rsi_FOD_matrix(Q,X,lADC,tADC)
%
%  Purpose: generate fiber orientation density forward matrix
%      
%  Required Input:
%    Q: [nm x 3] matrix of q values
%       diffusion directions, scaled by sqrt(b-value)
%    X: [nr x 3] matrix of FOD reconstruction points
%    lADC: longtitudinal ADC
%    tADC: transverse ADC 
%
%  Output:
%    R: [nm x nr] foward FOD signal matrix
% 
% Created:  05/08/12 by Don Hagler
% Last Mod: 10/27/12 by Don Hagler
%

% based on code by Nate White

if ~mmil_check_nargs(nargin,1), return; end;

nm = size(Q,1);
nr = size(X,1);
b = sum(abs(Q).^2,2);
U = zeros(nm,3);
for i = 1:nm
  U(i,:) = Q(i,:)/(sqrt(b(i))+eps);
end

R = zeros(nm,nr);
for i = 1:nr
  R(:,i) = exp(-b.*((lADC-tADC)*(U*X(i,:)').^2+tADC));
end

return;

