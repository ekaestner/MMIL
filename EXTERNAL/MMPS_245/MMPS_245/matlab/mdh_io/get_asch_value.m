function val = get_mdh_value(fid,tag)

fseek(fid,0,-1);

val = 0;
str = fgetl(fid);
while (~feof(fid))
  [tok,rest] = strtok(str);
  if (strcmp(tok,tag))
    [tok2,rest] = strtok(rest);
    [tok3,rest] = strtok(rest);
    if (tok(1)~='t')
      val = sscanf(tok3,'%f');
    else
      val = tok3;
    end
%    tok
%    tok2
%    tok3
%    val
    return;
  end
  str = fgetl(fid);
end
