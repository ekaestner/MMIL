function [fname_data,fname_grad,fname_bval]=...
  DTI_MMIL_Export_Data_for_FSL(ContainerPath,varargin)
%function [fname_data,fname_grad,fname_bval]=...
%  DTI_MMIL_Export_Data_for_FSL(ContainerPath,[options])
%
% Purpose:
%   Concatenates raw data files (specified by snums)
%     and writes diffusion gradient table as text file (corrected for motion) for use with FSL
%
% Usage:
%  DTI_MMIL_Export_Data(ContainerPath,'key1', value1,...);
%
% Required Input Parameters:
%   ContainerPath: full path of directory containing processed diffusion data
%                  (mgh/mgz format)
%
% Optional Parameters specifying which data to load:
%   'snums': list of scan numbers to concatenate and analyze
%     if empty (or unspecified), use all DTI scans in container
%     {default = []}
%   'infix': if empty, will look for files like 'DTI1.mgz'
%     otherwise, input file will be sprintf('DTI%d_%s.mgz',snum,infix)
%     example infix = 'gruw_reg_mc_ecc_B0uw_iso'
%     {default = []}
%   'revflag': [0|1|2] specify whether to use non-rev or rev data
%     if revflag=0, use non-rev data
%     if revflag=1, use rev data
%       rev scans have names like 'DTI1_rev.mgz'
%     if revflag=2, use concatenated non-rev and rev data
%     {default = 0}
%   'min_ndirs': require at least this many diffusion directions to be valid
%     {default = 6}
%   'min_bval': minimum b value a scan must have to be included
%     {default = 1}
%   'max_bval': maximum b-value used in tensor fit
%     {default = Inf}
%   'flex_flag': [0|1] DTI_flex scans included in tensor fit
%     {default = 0}
%   'min_nb0': minimum number of b=0 image to be included
%     {default = 1}
%   'nob0_flag': [0|1] toggle exclusion of b=0 images from fitting
%    if 1, multiple b-values are required
%    also, b=0 images are still used for between image scaling
%    {default = 0}
%
% Optional Parameters specifying output
%  'outdir': output directory
%    {default: ContainerPath}
%  'sort_b0vols_flag': [0|1] whether to sort b=0 images to beginning of volume
%    {default = 0}
%  'split_flag': [0|1] whether to export data as separate scans if there are two DTI scans 
%                 if 0, data is to be concatenated then 'sort_b0vols_flag' is set to 1
%    {default = 1}
%  'forceflag': [0|1] whether to run conversions even if output files exist
%    {default = 0}
%
% Optional Output:
%   fname_data: full path of output data file
%   fname_grad: full path of output gradient table file
%   fname_bval: full path of output bvalue file
%
% see also: DTI_MMIL_Export_DTIstudio
% see also: DTI_MMIL_Export_Data
%
% Created:  11/22/11 by Vijay Venkatraman
% Prev Mod: 11/21/13 by Don Hagler
% Last Mod: 06/23/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parse input parameters
if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'snums',[],[],...
  'infix',[],[],...
  'revflag',0,[0,1,2],...
  'min_nb0',1,[],...
  'min_ndirs',6,[],...
  'min_bval',1,[],...
  'max_bval',Inf,[100,Inf],...
  'flex_flag',false,[false true],...
  'nob0_flag',false,[false true],...
...
  'outdir','exportDTIforFSL',[],...
  'sort_b0vols_flag',false,[false true],...
  'smf',10^-5,[10^-100,10^-1],...
  'split_flag',true,[false true],...
  'forceflag',false,[false true],...  
});

if isempty(ContainerPath)
  fprintf('%s: ERROR: ContainerPath is empty\n',mfilename);
  errcode = 1;
  return;
elseif ~exist(ContainerPath,'file')
  fprintf('%s: ERROR: ContainerPath %s not found\n',mfilename,ContainerPath);
  errcode = 1;
  return;
end;

if mmil_isrelative(parms.outdir)
  parms.outdir = [ContainerPath '/' parms.outdir];
end;
mmil_mkdir(parms.outdir);

if isempty(parms.snums)
  tags = {'snums','revflag','min_bval','flex_flag','min_ndirs','min_nb0'};
  args = mmil_parms2args(parms,tags);
  [ScanInfo, SessInfo,errcode]= DTI_MMIL_Get_ScanInfo(ContainerPath,args{:});
  parms.snums = SessInfo.snums_DT;
end;

if ~parms.split_flag
  parms.sort_b0vols_flag = 1;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copy the DT measures to the output directory 
tags = {'snums','infix','revflag','min_bval','max_bval','flex_flag',...
  'min_ndirs','min_nb0','nob0_flag'};
args = mmil_parms2args(parms,tags);
[fstem,snums] = DTI_MMIL_Set_DT_fstem(ContainerPath,args{:});
if isempty(fstem)
  errcode = 1;
  return;
end;
fname_FA = sprintf('%s_FA.mgz',fstem);
fname_MD = sprintf('%s_MD.mgz',fstem);
fname_TD = sprintf('%s_TD.mgz',fstem);
fname_LD = sprintf('%s_LD.mgz',fstem);
fname_b0 = sprintf('%s_b0.mgz',fstem);
meas_outdir = sprintf('%s/DTmeas',parms.outdir);
if ~parms.forceflag
  all_exist = 1;
  outlist = {'nodif','FA','MD','RD','L1'};
  for i=1:length(outlist)
    fname_out = sprintf('%s/data_%s.nii.gz',meas_outdir,outlist{i});
    if ~exist(fname_out,'file')
      all_exist = 0;
    end;
  end;
  fname_out = sprintf('%s/nodif_brain_mask.nii.gz',meas_outdir);
  if ~exist(fname_out,'file')
    all_exist = 0;
  end;
else
  all_exist = 0;
end;
if ~all_exist
  cmd = sprintf('fsl_dti_measlist.csh %s %s %s %s %s %s %d',fname_FA,...
        fname_MD,fname_TD,fname_LD,fname_b0,meas_outdir,parms.forceflag);
  fprintf('%s: exporting DT measures to output directory for FSL...\n',mfilename);
  fprintf('%s\n',cmd);
  [status,r] = mmil_unix(cmd);
  if status
    fprintf('%s: cmd %s failed:\n%s\n',mfilename,cmd,r);
    errcode = 1;
    return;
  end;
end;

% copy the DT data to the output directory
if parms.split_flag
  for i= 1:length(parms.snums)
    s = parms.snums(i);
    snums_dir = sprintf('%s/DTI%d',parms.outdir,s);
    mmil_mkdir(snums_dir);  
    % assign the dti output directory
    parms.outdir_dti = snums_dir;
    % make fname_data a full path    
    parms.fname_data = sprintf('%s/data',snums_dir);
    fname_data{i} = parms.fname_data; 
    % make fname_grad a full path
    parms.fname_grad = [snums_dir '/bvecs.txt'];
    fname_grad{i} = parms.fname_grad;  
    % make fname_bval a full path
    parms.fname_bval = [snums_dir '/bvals.txt'];
    fname_bval{i} = parms.fname_bval;      
    % subfunction
    savedata(ContainerPath,s,parms);
  end;  
else
  snums_dir = parms.outdir;
  mmil_mkdir(snums_dir);  
  % assign the dti output directory
  parms.outdir_dti = snums_dir;  
  % make fname_data a full path
  parms.fname_data = sprintf('%s/data',snums_dir);
  fname_data = parms.fname_data;   
  % make fname_grad a full path
  parms.fname_grad = [snums_dir '/bvecs.txt'];
  fname_grad = parms.fname_grad;
  % make fname_bval a full path
  parms.fname_bval = [snums_dir '/bvals.txt'];
  fname_bval = parms.fname_bval;
  % subfunction 
  savedata(ContainerPath,parms.snums,parms);
end;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
function savedata(ContainerPath,s,parms)
  fname_out = sprintf('%s/data.nii.gz',parms.outdir_dti);
  if exist(parms.fname_bval,'file') &&...
     exist(parms.fname_grad,'file') &&...
     exist(fname_out,'file') &&...
     ~parms.forceflag
    return;
  end;

  % load concatenated volume and gradient directions
  [vol,M,qmat,bvals,TEs]=DTI_MMIL_Load_Data(ContainerPath,...
    'snums',s,'infix',parms.infix,'revflag',parms.revflag,...
    'min_nb0',parms.min_nb0,'min_ndirs',parms.min_ndirs,...
    'min_bval',parms.min_bval,...
    'flex_flag',parms.flex_flag);
  if isempty(vol), return; end;

  % set bval to 0 if length of qvec is zero
  qlength = sqrt(sum(qmat.^2,2));
  i_b0 = find(qlength<parms.smf);
  nb0 = length(i_b0);
  bvals(i_b0) = 0;

  % sort b=0 images to beginning
  if parms.sort_b0vols_flag
    [bvals,ind] = sort(bvals);
    qmat = qmat(ind,:);
    vol = vol(:,:,:,ind);
  end;

  % write bvals text file
  if ~exist(parms.fname_bval,'file') || parms.forceflag
    fid = fopen(parms.fname_bval,'wt');
    if fid==-1, error('failed to open file %s for writing',parms.fname_bval); end;
    fprintf(fid,'%d\n',bvals);
    fclose(fid);
  end;

  % write gradtable text file
  % reorient into FSL format LAS from RAS
  qmat(nb0+1:end,1) = -qmat(nb0+1:end,1);
  if ~exist(parms.fname_grad,'file') || parms.forceflag
    fid = fopen(parms.fname_grad,'wt');
    if fid==-1, error('failed to open file %s for writing',parms.fname_grad); end;
    for i=1:size(qmat,1)
      fprintf(fid,'%s\n',strtrim(sprintf('%0.6f ',qmat(i,:))));
    end;
    fclose(fid);
  end;

  % copy the DT data to the output directory
  fdata = sprintf('%s.mgz',parms.fname_data);
  fs_save_mgh(vol,fdata,M);

  % reorient to LAS orientation for FSL
  fdata_LAS = sprintf('%s_LAS.mgz',parms.fname_data);
  epi_reorient_vol(fdata,fdata_LAS,'LAS');
  cmd = sprintf('rm %s',fdata);
  [status,r] = mmil_unix(cmd);
  if status
    fprintf('%s: cmd %s failed:\n%s\n',mfilename,cmd,r);
    errcode = 1;
    return;
  end;
  fdata = sprintf('%s.mgz',parms.fname_data);  
  cmd = sprintf('mv %s %s',fdata_LAS,fdata);
  [status,r] = mmil_unix(cmd);
  if status
    fprintf('%s: cmd %s failed:\n%s\n',mfilename,cmd,r);
    errcode = 1;
    return;
  end;

  % convert MGH format to NIFTIGZ format
  cmd = sprintf('fsl_dti_data.csh %s %s %d',...
    fdata,parms.outdir_dti,parms.forceflag);
  fprintf('%s: exporting DT data to output directory for FSL...\n',mfilename);
  fprintf('%s\n',cmd);
  [status,r] = mmil_unix(cmd);
  if status
    fprintf('%s: cmd %s failed:\n%s\n',mfilename,cmd,r);
    errcode = 1;
    return;
  end;
return;


