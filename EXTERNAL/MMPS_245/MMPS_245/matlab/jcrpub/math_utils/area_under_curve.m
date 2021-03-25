function [a] = area_under_curve(x,y)

if length(x) ~= length(y)
   error('X and Y must have the same length.');
end

if length(x) <= 1,
   error('Length of X and Y must be greater than 1');
end

a = 0;
for ii=2:length(x)
   yval = 0.5*(y(ii)+y(ii-1));
   a = a + yval*abs(x(ii)-x(ii-1));
end
