function y = mmil_gauss2_eval(x1,x2,beta)
%function y = mmil_gauss2_eval(x1,x2,beta)
%
% Purpose: evaluate 2D Gaussian function of the form
%   Y = a*exp(-(p*(x1-u1).^2 + 2*q*(x1-u1).*(x2-u2) + r*(x2-u2).^2))
%   with p = cos(t)^2/(2*s1^2) + sin(t)^2/(2*s2^2)
%        q = -sin(2*t)/(4*s1^2) + sin(2*t)/(4*s2^2)
%        r = sin(t)^2/(2*s1^2) + cos(t)^2/(2*s2^2)
%
% Required Input:
%   x1: vector or matrix of values for independent variable 1
%   x2: vector or matrix of values for independent variable 2
%     size of x2 must match x1
%   beta: vector of Gaussian fit parameters [a,u1,s1,u2,s2,t]
%
% Output:
%   y: output vector or matrix with size of x1
%
% Created:  08/14/12 by Don Hagler
% Last Mod: 08/14/12 by Don Hagler
%

Y = [];
if ~mmil_check_nargs(nargin,3), return; end;

if length(size(x1)) ~= length(size(x2)) || ...
   any(size(x1)~=size(x2))
  error('size of x1 does not match x2');
end;

% get fit parameters from beta vector
if numel(beta)>=1
  a = beta(1);
else
  a = 1;
end;
if numel(beta)>=3
  u1 = beta(2);
  s1 = beta(3);
else
  u1 = 0;
  s1 = 1;
end;
if numel(beta)>=5
  u2 = beta(4);
  s2 = beta(5);
else
  u2 = u1;
  s2 = s1;
end
if numel(beta)>5
  t = beta(6);
else
  t = 0;
end;

% calculate function coefficients from fit parameters
p = cos(t)^2/2/s1^2 + sin(t)^2/2/s2^2;
q = -sin(2*t)/4/s1^2 + sin(2*t)/4/s2^2;
r = sin(t)^2/2/s1^2 + cos(t)^2/2/s2^2;

% evaluate 2D Gaussian function
y = a*exp(-(p*(x1-u1).^2 + 2*q*(x1-u1).*(x2-u2) + r*(x2-u2).^2));

return;

