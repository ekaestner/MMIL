function randval = rand_float(valrange)
%function randval = rand_float([minval,maxval])
%
%
%

if length(valrange)<2
  valrange(2) = valrange;
  valrange(1) = 0;
end;

randval = rand(1)*(valrange(2)-valrange(1)) + valrange(1);

