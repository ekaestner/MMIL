function rc_write_cond_info(cond_info,fname_conds,forceflag)
%function rc_write_cond_info(cond_info,fname_conds,forceflag)
%
% Created:  02/05/09 by Don Hagler
% Last Mod: 01/23/13 by Don Hagler
%

if ~mmil_check_nargs(nargin,2), return; end;
if ~exist('forceflag','var') || isempty(forceflag), forceflag=0; end;

if ~exist(fname_conds,'file') || forceflag
  fid = fopen(fname_conds,'wt');
  if fid==-1, error('failed to open file %s for writing',fname_conds); end;

  flist = fieldnames(cond_info);
  for f=1:length(flist)
    if f>1, fprintf(fid,','); end;
    fprintf(fid,'"%s"',flist{f});
  end;
  fprintf(fid,'\n');

  for i=1:length(cond_info)
    for f=1:length(flist)
      if f>1, fprintf(fid,','); end;
      fieldname = flist{f};
      val = cond_info(i).(fieldname);
      if ~isempty(findstr('offset',fieldname))
        if isempty(val)
          fprintf(fid,'0');
        else
          fprintf(fid,'%0.7f',val);
        end;
      elseif isfloat(val) & mmil_isint(val)
        fprintf(fid,'%d',val);
      elseif isfloat(val)
        fprintf(fid,'%0.7f',val);
      elseif ischar(val)
        fprintf(fid,'"%s"',val);            
      end;
    end;        
    fprintf(fid,'\n');
  end;
  fclose(fid);
end;

