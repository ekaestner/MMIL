function MMIL_Make_BEM_Surfs(ContainerPath,FSContainerPath,forceflag)
%function MMIL_Make_BEM_Surfs(ContainerPath,FSContainerPath,forceflag)
%
% Created:  06/14/09 by Don Hagler
% Last Mod: 06/16/16 by Don Hagler
%


if ~mmil_check_nargs(nargin, 1), return; end;
if ~exist('FSContainerPath','var') FSContainerPath=[]; end;
if ~exist('forceflag','var') | isempty(forceflag), forceflag=0; end;
if FSContainerPath==1 | FSContainerPath==0
  forceflag=FSContainerPath;
  FSContainerPath = [];
end;
ext = '.mgz';

fprintf('%s(''%s'',''%s'',%d)\n',mfilename,...
  ContainerPath,FSContainerPath,forceflag);

%fprintf('%s: making BEM surfs...\n',mfilename);

% BEM from T1
fname_in = sprintf('%s/MPR_res%s',ContainerPath,ext);
if exist(fname_in,'file')
  outdir = sprintf('%s/bem_T1',ContainerPath);
  mmil_makeBEMsurfs_fromT1(fname_in,outdir,forceflag);
else
  fname_in = sprintf('%s/hiFA_res%s',ContainerPath,ext);
  if exist(fname_in,'file')
    outdir = sprintf('%s/bem_T1',ContainerPath);
    mmil_makeBEMsurfs_fromT1(fname_in,outdir,forceflag);
  end;
end;

% BEM from aseg
if exist(fname_in,'file') & ~isempty(FSContainerPath)
  fname_aseg = sprintf('%s/mri/aparc+aseg.mgz',FSContainerPath);
  if ~exist(fname_aseg,'file')
    fprintf('%s: WARNING: file %s not found\n',mfilename,fname_aseg);
  else
    outdir = sprintf('%s/bem_aseg',ContainerPath);
    mmil_makeBEMsurfs_from_aseg(fname_in,fname_aseg,outdir,forceflag);
  end;
end;

% BEM from PD
fname_in = sprintf('%s/loFA_res%s',ContainerPath,ext);
if exist(fname_in,'file')
  outdir = sprintf('%s/bem_PD',ContainerPath);
  mmil_makeBEMsurfs_fromPD(fname_in,outdir,forceflag);

  % BEM from PD plus aseg with NFT
  if ~isempty(FSContainerPath)
    mmil_makeBEMsurfs_withNFT(ContainerPath,FSContainerPath,...
      'fstem_aseg','aparc+aseg',...
      'fname_brain',fname_in,...
      'input_bem_dir','bem_PD',...
      'output_bem_dir','bem_PD_NFT',...
      'outstem','PD_NFT',...
      'forceflag',forceflag);
  end;
end;

