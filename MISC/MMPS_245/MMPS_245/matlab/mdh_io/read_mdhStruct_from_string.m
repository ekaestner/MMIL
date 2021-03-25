function mdhStruct = read_mdhStruct_from_string(hdrstr)
%% Read mdh.asc file and return mdh fields

mdhStruct = struct();
matindx = 0;
endlnvec = strfind(hdrstr,char(10));
endlnvec = [0 endlnvec length(hdrstr)+1];
for linenum = 1:(length(endlnvec)-1)
  str = hdrstr((endlnvec(linenum)+1):(endlnvec(linenum+1)-1));
  eqpos = find(str=='=');
  if (length(str)>0)&(str(1)=='#')&(length(eqpos)>=1)
    hashpos = find(str=='#');
    rest = str((hashpos(end)+2):end);
    while ~isempty(rest)
      [lhs,rest] = strtok(rest);
      [junk,rest] = strtok(rest);
      [rhs,rest] = strtok(rest);
      if ~isempty(lhs) & ~isempty(rhs)
        lhs = subst_lhs_chars(lhs);
        if ~isnumeric(rhs)
          tmp = str2double(rhs);
          if ~isnan(tmp) rhs = tmp; end
        end
        mdhStruct = setfield(mdhStruct,lhs,rhs);
        if strmatch('adRM',lhs)
          if matindx==0 mdhStruct.adRM = zeros(3,3); end
          matindx = matindx+1;
          mdhStruct.adRM(matindx) = rhs;
        end
      end
    end
  elseif (length(str)>0)&(str(1)~='#')&(length(eqpos)==1)
    lhs = str(1:(eqpos-1));
    rhs = str((eqpos+2):end);
    quotepos = find(rhs=='"');
    if length(quotepos)>0
      rhs = rhs((quotepos(1)+1):(quotepos(end)-1));
    end
    lsqbpos = find(rhs=='[');
    rsqbpos = find(rhs==']');
    if (length(lsqbpos)==1)&(length(rsqbpos)==1)
      numvals = 0;
      rhsval = [];
      rest = rhs((lsqbpos+1):(rsqbpos-1));
      while ~isempty(rest)
        [tok,rest] = strtok(rest);
        if ~isempty(tok)
          numvals = numvals+1;
          rhsval = [rhsval tok ' '];
        end
      end
    else
      rhsval = rhs;
    end
    lhs = subst_lhs_chars(lhs);
    if ~isnumeric(rhsval)
      rhsval = deblank(rhsval);
      tmp = str2double(rhsval);
      if isnan(tmp) tmp =  str2num(rhsval); end
      if ~isempty(tmp) rhsval = tmp; end
    end
    mdhStruct = setfield(mdhStruct,lhs,rhsval);
  end
end
if matindx>0 mdhStruct(1).adRM = mdhStruct(1).adRM.'; end
