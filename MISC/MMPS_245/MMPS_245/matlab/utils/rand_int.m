function randval = rand_int(valrange)
%function randval = rand_int([minval,maxval])
%
%
%

if length(valrange)<2
  valrange(2) = valrange;
  valrange(1) = 0;
end;

randval = round(rand(1)*(valrange(2)-valrange(1)) + valrange(1));

