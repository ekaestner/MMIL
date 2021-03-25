function errcode = MMIL_Copy_BEM_Surfs(ContainerPath,FSContainerPath,BEMtype,forceflag)
%function errcode = MMIL_Copy_BEM_Surfs(ContainerPath,FSContainerPath,BEMtype,forceflag)
%
% Created:  05/13/10 by Don Hagler
% Last Mod: 06/16/14 by Don Hagler
%

errcode = 0;
if (~mmil_check_nargs(nargin, 2)), return; end;
if ~exist('BEMtype','var') | isempty(BEMtype), BEMtype='aseg'; end;
if ~exist('forceflag','var') | isempty(forceflag), forceflag=0; end;

fprintf('%s(''%s'',''%s'',''%s'',%d)\n',mfilename,...
  ContainerPath,FSContainerPath,BEMtype,forceflag);

bem_list = {'PD','PD_NFT','T1','aseg'};

if ~ismember(BEMtype,bem_list)
  error('invalid BEMtype: %s',BEMtype);
end;

PD_NFT_suffix = '_n10000_d2000'; %% todo: expose this as an option?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(ContainerPath,'dir')
  error('processed container %s not found',ContainerPath);
end;
if ~exist(FSContainerPath,'dir')
  error('FreeSurfer container %s not found',FSContainerPath);
end;

k = 0;
for b=1:length(bem_list)
  in_bem_dir = [ContainerPath '/bem_' bem_list{b}];
  if exist(in_bem_dir,'dir')
    k = k + 1;
    out_bem_dir = [FSContainerPath '/bem_' bem_list{b}];
    if exist(out_bem_dir,'dir') & ~isempty(dir([out_bem_dir '/*.tri']))
      if ~forceflag
        fprintf('%s: WARNING: output dir %s already exists, not replacing\n',...
          mfilename,out_bem_dir);
        continue;
      else
        fprintf('%s: WARNING: output dir %s already exists, deleting and replacing\n',...
          mfilename,out_bem_dir);    
        errcode = unix_cmd(['rm -r ' out_bem_dir]);
        if errcode, return; end;
      end;
    end;
    errcode=unix_cmd(['cp -rpd ' in_bem_dir ' ' out_bem_dir]);
    if errcode, return; end;
  end;
end;
if k==0
  fprintf('%s: WARNING: no bem directories found in %s\n',mfilename,ContainerPath);
  return;
else
  fprintf('%s: %d bem directories found in %s\n',mfilename,k,ContainerPath);
end;


out_bem_dir = [FSContainerPath '/bem'];
if exist(out_bem_dir,'dir') & ~isempty(dir([out_bem_dir '/*.tri']))
  if ~forceflag
    fprintf('%s: WARNING: output dir %s already exists, not replacing\n',...
      mfilename,out_bem_dir);
    return;
  else
    fprintf('%s: WARNING: output dir %s already exists, deleting and replacing\n',...
      mfilename,out_bem_dir);    
    errcode = unix_cmd(['rm -r ' out_bem_dir]);
    if errcode, return; end;
  end;
end;
[succ,msg] = mkdir(out_bem_dir);
if ~succ
  fprintf('%s: ERROR: failed to create output dir %s\n',mfilename,outdir);
  errcode = 1; return;
end;


if strcmp(BEMtype,'PD_NFT')
  in_bem_dir = [FSContainerPath '/bem_PD_NFT'];
  if exist(in_bem_dir,'dir')
    errcode = unix_cmd(sprintf('ln -s %s/brain%s.tri %s/brain.tri',in_bem_dir,PD_NFT_suffix,out_bem_dir));
    if errcode, return; end;
    errcode = unix_cmd(sprintf('ln -s %s/inner_skull%s.tri %s/inner_skull.tri',in_bem_dir,PD_NFT_suffix,out_bem_dir));
    if errcode, return; end;
    errcode = unix_cmd(sprintf('ln -s %s/outer_skull%s.tri %s/outer_skull.tri',in_bem_dir,PD_NFT_suffix,out_bem_dir));
    if errcode, return; end;
    errcode = unix_cmd(sprintf('ln -s %s/outer_scalp%s.tri %s/outer_scalp.tri',in_bem_dir,PD_NFT_suffix,out_bem_dir));
    if errcode, return; end;
  else
    fprintf('%s: WARNING: bem_PD_NFT dir %s not found\n',mfilename,in_bem_dir);
  end;
elseif strcmp(BEMtype,'PD')
  in_bem_dir = [FSContainerPath '/bem_PD'];
  if exist(in_bem_dir,'dir')
    errcode = unix_cmd(sprintf('ln -s %s/inner_skull4.tri %s/inner_skull.tri',in_bem_dir,out_bem_dir));
    if errcode, return; end;
    errcode = unix_cmd(sprintf('ln -s %s/outer_skull4.tri %s/outer_skull.tri',in_bem_dir,out_bem_dir));
    if errcode, return; end;
    errcode = unix_cmd(sprintf('ln -s %s/outer_scalp4.tri %s/outer_scalp.tri',in_bem_dir,out_bem_dir));
    if errcode, return; end;
  else
    fprintf('%s: WARNING: bem_PD dir %s not found\n',mfilename,in_bem_dir);
  end;
else
  in_bem_dir = [FSContainerPath '/bem_T1'];
  if exist(in_bem_dir,'dir')
    errcode = unix_cmd(sprintf('ln -s %s/outer_skull4.tri %s/outer_skull.tri',in_bem_dir,out_bem_dir));
    if errcode, return; end;
    errcode = unix_cmd(sprintf('ln -s %s/outer_scalp4.tri %s/outer_scalp.tri',in_bem_dir,out_bem_dir));
    if errcode, return; end;
  else
    fprintf('%s: WARNING: bem_T1 dir %s not found\n',mfilename,in_bem_dir);
  end;
  if strcmp(BEMtype,'aseg')
    in_bem_aseg_dir = [FSContainerPath '/bem_aseg'];
    if exist(in_bem_aseg_dir,'dir')
      errcode = unix_cmd(sprintf('ln -s %s/inner_skull_aseg.tri %s/inner_skull.tri',...
        in_bem_aseg_dir,out_bem_dir));
      if errcode, return; end;
    else
      fprintf('%s: WARNING: bem_aseg dir %s not found\n',mfilename,in_bem_aseg_dir);
    end;
  else
    if exist(in_bem_dir,'dir')
      errcode = unix_cmd(sprintf('ln -s %s/inner_skull4.tri %s/inner_skull.tri',in_bem_dir,out_bem_dir));
      if errcode, return; end;
    end;
  end;
end;

return;



function errcode = unix_cmd(cmd)
  errcode = 0;
  [s,r] = unix(cmd);
  if s
    fprintf('%s: ERROR: cmd %s failed:\n%s',mfilename,cmd,r);
    errcode = 1;
  end;
return;
