function  D = StripQuotes(Din)

D = Din;
for i = 1:size(D,1)
  for j = 1:size(D,2)
    if ischar(D{i,j})
      s = D{i,j};
      if s(1) == '"' & s(end) == '"'
        D{i,j} = s(2:end-1);
      end
    end
  end
end
