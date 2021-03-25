function txtnum = ConvertXLS2mat(xlsfname,matfname)

[numeric,txt] = xlsread(xlsfname);
xls_txtnum = {};
xls_num = zeros(size(txt));
% Fix stupid xlsread problem with  text vs numerical data
for i = 1:size(txt,1)
  if mod(i,10)==0 fprintf(1,'Row %d of %d\n',i,size(txt,1)); end
  for j = 1:size(txt,2)
    if (j-1)>=1 & (j-1)<=size(numeric,2) & isfinite(numeric(i,j-1))
      f = numeric(i,j-1);
    else 
      f = str2num(txt{i,j});
    end
    if ~isempty(f)
      xls_txtnum{i,j} = f;
      xls_num(i,j) = f;
    else
      xls_txtnum{i,j} = txt{i,j};
      xls_num(i,j) = NaN;
    end
  end
end

save(matfname,'xls_txtnum','xls_num');

