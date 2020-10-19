function errcode = MMIL_Register_and_Resample_MRI_Volumes(ContainerPath,varargin);
%function errcode = MMIL_Register_and_Resample_MRI_Volumes(ContainerPath,[options]);
%
% Required Input:
%   ContainerPath: full path of MRIPROC Container
%
% Optional Input:
%  'atlasflag': [0|1] whether to resample output in atlas space
%              (rigid-body registration only)
%    {default = 0}
%  'nativeflag': [0|1] if atlasflag=0, whether to keep native resolution
%     otherwise, resample to 1mm, 256^3, LIA
%    {default = 0}
%  'rawQCflag': [0|1] whether to require manual raw QC
%     Looks for QC file in /home/<USER>/ProjInfo/<ProjID>/<ProjID>_RawQC.mat
%    {default = 0}
%  'ProjID': Project Name (required if rawQCflag=1)
%    {default = []}
%  'minmax': vector of minimum and maximum numbers of scans per session to register
%     If empty, use all available scans
%    {default = []}
%  'scantypes': cell array of scan types to register
%    {default = {'MPR','XetaT2','FLASHhi','FLASHlo','MEDIChi','MEDIClo'}}
%  'refT1flag': which type of T1 series ('MPR' or 'hiFA') to use as reference
%     i.e. which one gets registered to atlas and to which other scans are registered
%     0=MPR; 1=hiFA; 2=Either (prefer MPR); 3=Either (prefer hiFA)
%    {default = 2}
%  'gradunwarp_flag': [0|1|2] whether to correct for gradient warping
%     if 2, processing is aborted if gradwarp info is missing
%    {default = 1}
%  'wmbc_flag': [0|1] whether images were corrected for intensity bias-field
%     using white matter segmentation and sparse smoothing
%    {default = 0}
%  'image_range': range of values to clip output image
%    only applies if wmbc_flag = 1
%    {default = [0,255]}
%  'T2w_image_range': range of values to clip output for T2-weighted image
%    {default = [0,2500]}
%  'nu_flag': [0|1] whether images were corrected for intensity bias-field 
%     using nu (N3) correction
%     NOTE: if wmbc_flag = 1, nu_flag is irrelevant
%    {default = 0}
%  'atlasdir': full path of atlas directory
%    {default =  [getenv('MMPS_DIR') '/atlases']}}
%  'atlasname': name of atlas file (omit .mat extension)
%     full path or relative to atlasdir
%    {default = 'T1_Atlas/T1_atlas'}
%  'oldreg_flag': [0|1] use old registration mriRegister
%     otherwise use mmil_reg in mmil_average_volumes
%    {default = 0}
%  'prereg_flag': [0|1] whether to use reg for rigid registration
%     before running dct registration to atlas
%    {default = 1}  
%  'mask_smooth1': smoothing sigma (voxels) for initial fill step
%    {default: 20}
%  'mask_thresh1': threshold applied to mask after first smoothing step
%    {default: 0.5}
%  'mask_smooth2': smoothing sigma (voxels) for second dilation step
%    {default: 40}
%  'mask_thresh2': threshold applied to mask after second smoothing step
%    {default: 0.2}
%  'mask_smooth3': smoothing sigma (voxels) for third dilation step
%    {default: 10}
%  'ext': file extension (i.e. '.mgh' or '.mgz')
%    {default = '.mgz'}
%  'forceflag': [0|1] whether to overwrite existing output
%    {default = 0}
% 
% Created:  03/18/08 by Don Hagler
% Prev Mod: 07/14/17 by Don Hagler
% Last Mod: 08/06/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% todo: for MEDIChi and MEDIClo, they are preregistered (same scan)
%%       so only need to register MEDIChi to MPR and then apply that M to MEDIClo
%%       except subsequent scans need to be registered to the new ref

%% todo: for MEDIChi and MEDIClo, apply registration for each echo time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errcode = 0;
if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'atlasflag',false,[false true],...
  'nativeflag',false,[false true],...
  'rawQCflag',false,[false true],...
  'ProjID',[],[],...
  'minmax',[1 Inf],[1 Inf],...
  'scantypes',{'MPR','XetaT2','FLASHhi','FLASHlo','MEDIChi','MEDIClo'},[],...
  'refT1flag',2,[0:3],...
  'wmbc_flag',false,[false true],...
  'image_range',[0,255],[],...
  'T2w_image_range',[0,2500],[],...
  'nu_flag',false,[false true],...
  'gradunwarp_flag',1,[0:2],...  
  'ext','.mgz',{'.mgh','.mgz'},...
  'forceflag',false,[false true],...
...
  'cleanupflag',true,[false true],...
  'hatl_offset_RAS',[0 -15 -10],[],... % displace resampled image from atlas (preserve nose)
  'mask_smooth',15,[],...
  'atlasdir',[],[],...
  'atlasname','T1_Atlas/T1_atlas',[],...
  'oldreg_flag',false,[false true],...
  'prereg_flag',true,[false true],...
  'mask_smooth1', 20,[0 100],...
  'mask_thresh1', 0.5,[0 1],...
  'mask_smooth2', 40,[0 100],...
  'mask_thresh2', 0.2,[0 1],... 
  'mask_smooth3', 10,[0 100],...  
  'nu_scantypes',{'MPR','XetaT2','FLASHhi','MEDIChi'},[],...
  'wmbc_scantypes',{'MPR','XetaT2','FLASHhi','MEDIChi'},[],...
  'T2w_scantypes',{'XetaT2'},[],...
...
  'ContainerPath',ContainerPath,[],...
});

if isempty(parms.atlasdir)
  parms.atlasdir = [getenv('MMPS_DIR') '/atlases'];
end;

parms.fname_atl = sprintf('%s/%s.mgz',parms.atlasdir,parms.atlasname);
parms.fname_atl_mask = sprintf('%s/%s_mask.mgz',parms.atlasdir,parms.atlasname);
parms.fname_atl_maskBroad = sprintf('%s/%s_maskBroad.mgz',parms.atlasdir,parms.atlasname);

% check fname_atl
if ~exist(parms.fname_atl,'file')
  error('%s: %s file missing',mfilename,parms.fname_atl);
end;

% check fname_atl_mask
if ~exist(parms.fname_atl_mask,'file')
   error('%s: %s file missing',mfilename,parms.fname_atl_mask);
end;

% check fname_atl_maskBroad
if ~exist(parms.fname_atl_maskBroad,'file')
  error('%s: %s file missing',mfilename,parms.fname_atl_maskBroad);
end;

parms.scantypes = unique(parms.scantypes);

if ~parms.wmbc_flag
  parms.image_range = [];
end;

M_hatl_offset = mmil_construct_M('trans',parms.hatl_offset_RAS);

if ~exist(ContainerPath,'dir')
  error('ContainerPath %s not found',ContainerPath);
end;
ScanInfo = [];
M_atl_to_subj = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ContainerInfo,errcode] = MMIL_Load_ContainerInfo(ContainerPath);
if errcode~=0, return; end;
ScanInfo = ContainerInfo.ScanInfo;

% load QC information
if parms.rawQCflag
  MPRflag = 0; FLASHhiflag = 0; % which series raw QC is done on
  qcfile = sprintf('%s/ProjInfo/%s/%s_RawQC.mat',...
    getenv('HOME'),parms.ProjID,parms.ProjID);
  if ~exist(qcfile,'file')
    fprintf('ERROR: RAWQC file %s not found\n',qcfile);
    errcode = 1;
    return;
  end
  load(qcfile);
  SubjIDs = {qcinfo.SubjID}; StudyDates = cell2mat({qcinfo.studydate});
  VisitIDs = {qcinfo.VisitID};
  [tmp_path,ContainerDir,tmp_ext] = fileparts(ContainerPath);
  ContainerDir = [ContainerDir tmp_ext];
  VisitID = char(regexp(ContainerDir,'(?<=MRIPROC_).+(?=_\d{8})','match'));
  StudyDate = str2num(char(regexp(ContainerDir,'(?<=_)\d{8}(?=\.)','match')));
  S_ind = find(strcmp(VisitIDs,VisitID) & StudyDates==StudyDate);
  if isempty(S_ind)
    fprintf('%s: ERROR: VisitID %s, StudyDate %d not found in QC file\n', ...
      mfilename,VisitID,StudyDate);
    errcode = 1;
    return;
  end
  SubjID = qcinfo(S_ind(1)).SubjID;
  qcscans = {};
  min_s = MinMax(1); max_s = MinMax(2);
  for s=1:length(S_ind) % for each series QC'd (3 step QC)
      qcscans(end+1).qc = qcinfo(S_ind(s)).qc;
      qcscans(end).fname = char(qcinfo(S_ind(s)).input);
  end
  MPRseries = regexp({qcscans.fname},'MPR'); MPRseries = [MPRseries{:}];
  hiFASeries = regexp({qcscans.fname},'FLASHhi'); hiFASeries = [hiFASeries{:}];
  switch parms.refT1flag
    case 0 % require MPR as ref
      if ~isempty(MPRseries), MPRflag = 1; end;
    case 1 % require hiFA as ref
      if ~isempty(hiFASeries), FLASHhiflag = 1; end;
    case 2 % prefer MPR as ref, allow hiFA
      if ~isempty(MPRseries)
          MPRflag = 1;
      elseif ~isempty(hiFASeries)
          FLASHhiflag = 1;
      end
    case 3 % prefer hiFA as ref, allow MPR
      if ~isempty(hiFASeries)
          FLASHhiflag = 1;
      elseif ~isempty(MPRseries)
          MPRflag = 1;
      end
  end;
  allqc = cell2mat({qcscans.qc});
  S_ind = find(allqc == 1 | allqc ==2); % omit QC = 0 or 3
  if isempty(S_ind) | length(S_ind) < min_s %if all bad or not qc'd
    fprintf('%s: WARNING: Not enough series w/ OK QC for Subject %s, StudyDate %d\n',...
      mfilename,SubjID,StudyDate);
    return;
  end
  qcscans = qcscans(S_ind);
  if length(S_ind) > max_s    % if more series than max, take max # of series w/ best QC 
    allqc = cell2mat({qcscans.qc});
    [b ix] = sort(allqc);
    qcscans = qcscans(ix);
    qcscans = qcscans(1:max_s);
  end
end

% look for files of each scantype
nscans = zeros(length(parms.scantypes),1);
for t=1:length(parms.scantypes)
  fnamestem = parms.scantypes{t};
  tmpinfo = mmil_getfield(ScanInfo,fnamestem);
  if isempty(tmpinfo)
    continue;
  else
    if parms.rawQCflag  % uses rawQC info, minmax
      n = regexp({qcscans.fname},fnamestem); 
      nscans(t) = sum(cell2mat(n));
    else
      nscans(t) = length(find([tmpinfo.valid]==1));
    end;
  end;
end;

% look for reference scan
ref_scantype = [];
switch parms.refT1flag
  case 0 % require MPR as ref
    allowed_scantypes = {'MPR'};
  case 1 % require hiFA as ref
    allowed_scantypes = {'FLASHhi','MEDIChi'};
  case 2 % prefer MPR as ref, allow hiFA
    allowed_scantypes = {'MPR','FLASHhi','MEDIChi'};
  case 3 % prefer hiFA as ref, allow MPR
    allowed_scantypes = {'FLASHhi','MEDIChi','MPR'};
end;
ind = find(ismember(upper(parms.scantypes),upper(allowed_scantypes)));
if isempty(ind)
  error('refT1flag = %d, but %s not included in scantypes',parms.refT1flag,...
    sprintf('''%s'' ',allowed_scantypes));
end;
for a=1:length(allowed_scantypes)
  ind = find(strcmp(upper(parms.scantypes),upper(allowed_scantypes{a})));
  if isempty(ind) || nscans(ind)<1
    continue;
  else
    ref_scantype = allowed_scantypes{a};
    break;
  end;
end;
if isempty(ref_scantype)
  fprintf('%s: ERROR: missing reference scan\n',mfilename);
  errcode = 1;
  return;
end;

% register to atlas and create brain mask for reference
fname_ref = []; fname_regmat = []; fname_mask = []; suffix = 'none';
if parms.atlasflag || sum(nscans)>1
  fnamestem = ref_scantype;
  switch fnamestem
    case {'FLASHhi','MEDIChi'}
      outstem_mask = 'hiFA';
    case {'FLASHlo','MEDIClo'}
      outstem_mask = 'loFA';
    case {'XetaT2'}
      outstem_mask = 'T2w';
    otherwise %% todo: set MPR to T1w?
      outstem_mask = fnamestem;
  end; 

  % determine correct file name suffix
  suffix = set_suffix(parms,ScanInfo,fnamestem);
  if strcmp(suffix,'none')
    fprintf('%s: ERROR: no reference scan found\n',mfilename);
    errcode = 1;
    return;
  end;

  % determine which file is to be the reference
  tmpinfo = mmil_getfield(ScanInfo,fnamestem);
  for i=1:length(tmpinfo)
    fname = sprintf('%s/%s%d%s%s',ContainerPath,fnamestem,i,suffix,parms.ext);
    [fpath,fstem,fext] = fileparts(fname);
    if parms.rawQCflag
      n = regexp({qcscans.fname},sprintf('%s%d',fnamestem,i)); n = [n{:}];
      if ((MPRflag && strcmp(fnamestem,'MPR')) ||...
          (FLASHhiflag && strcmp(fnamestem,'FLASHhi'))) &&...
          isempty(n) %
%        fprintf('%s: BAD QC: skipping file %s%s\n',mfilename,fstem,fext);
        continue;
      end
    end
    if exist(fname,'file')
      fname_ref = fname;
      break;
    end;
  end;
  if isempty(fname_ref)
    fprintf('%s: ERROR: no valid reference scans found\n',mfilename);
    errcode = 1;
    return;
  end;
  
  % register to atlas and create brain mask
  [fpath,fstem] = fileparts(fname_ref);
  fname_regmat = sprintf('%s/%s_reg2atl.mat',ContainerPath,fnamestem);
  fname_mask = sprintf('%s/%s_mask%s',ContainerPath,fstem,parms.ext);
  fname_ref_res = sprintf('%s/%s_res%s',ContainerPath,fstem,parms.ext); 
  fname_mask_res = sprintf('%s/%s_res_mask%s',ContainerPath,outstem_mask,parms.ext);
  if ~exist(fname_ref_res,'file') || ~exist(fname_mask_res,'file') || parms.forceflag
    if ~exist(fname_mask,'file') || ~exist(fname_regmat,'file') ||...
       (~exist(fname_mask_res,'file') && parms.prereg_flag) ||...
       parms.forceflag
      if parms.prereg_flag
        % initial rigid registration to atlas (output is LIA,256^3 and 1mm^3) and this includes the offset
        fprintf('%s: rigid pre-registration to atlas of %s...\n',...
          mfilename,fname_ref);
        fname_ref_prereg = sprintf('%s/%s_prereg%s',...
          ContainerPath,fstem,parms.ext); 
        M_atl_to_subj_prereg = initial_rigidreg(fname_ref,fname_ref_prereg,...
          M_hatl_offset,parms);
        % registration to atlas to create brain mask
        fprintf('%s: nonlinear registration to atlas of %s...\n',...
          mfilename,fname_ref_prereg);
        fname_mask_prereg = sprintf('%s/%s_prereg_mask%s',...
          ContainerPath,outstem_mask,parms.ext);
        [volmask_prereg,M_atl_to_subj_postreg] = mmil_dct_brainmask([],...
          'fname_in',fname_ref_prereg,...
          'fname_mask',fname_mask_prereg,...
          'atlasdir',parms.atlasdir,'atlasname',parms.atlasname,...
          'smooth',parms.mask_smooth,'forceflag',1);
        [volmask_prereg,M_prereg] = ctx_ctx2mgh(volmask_prereg);
        % combine M_atl_to_subj_prereg and M_atl_to_subj_postreg
        M_atl_to_subj = M_atl_to_subj_prereg * M_atl_to_subj_postreg;
        % resample ref to res
        fprintf('%s: resampling ref volume...\n',mfilename);
        [M_res,volsz_res] = mmil_load_mgh_info(fname_ref_prereg,parms.forceflag);
        [vol_ref,M_ref,mr_parms,volsz_ref] = fs_load_mgh(fname_ref);
        vol_ref_res = mmil_resample_vol(vol_ref,M_ref,...
          'M_ref',M_res,'nvox_ref',volsz_res(1:3),...
          'M_reg',M_atl_to_subj);
        fs_save_mgh(vol_ref_res,fname_ref_res,M_res);
        % resample volmask_prereg to volmask_res
        fprintf('%s: resampling mask volumes...\n',mfilename);
        volmask_res = mmil_resample_vol(volmask_prereg,M_prereg,...
          'M_ref',M_res,'nvox_ref',volsz_res(1:3),...
          'M_reg',M_atl_to_subj_postreg);
        volmask_res = 1.0*(volmask_res>0.5);
        fs_save_mgh(volmask_res,fname_mask_res,M_res);
        % resample mask back to reference volume
        volmask_ref = mmil_resample_vol(volmask_res,M_res,...
          'M_ref',M_ref,'nvox_ref',volsz_ref(1:3),...
          'M_reg',inv(M_atl_to_subj));
        volmask_ref = 1.0*(volmask_ref>0.5);
        fs_save_mgh(volmask_ref,fname_mask,M_ref);
      else
        fprintf('%s: nonlinear registration to atlas of %s...\n',...
          mfilename,fname);
        [volmask_ref,M_atl_to_subj] = mmil_dct_brainmask([],...
          'fname_in',fname_ref,...
          'fname_mask',fname_mask,...
          'atlasdir',parms.atlasdir,'atlasname',parms.atlasname,...
          'smooth',parms.mask_smooth,'forceflag',1);
        % overwriting the existing file after threshold
        [volmask_ref,M_ref] = ctx_ctx2mgh(volmask_ref);
        volmask_ref = 1.0*(volmask_ref>0.5);
        fs_save_mgh(volmask_ref,fname_mask,M_ref);
        M_atl_to_subj = M_atl_to_subj*M_hatl_offset;  
      end;
      clear volmask volmask_res volmask_ref;    
      save(fname_regmat,'M_atl_to_subj','fname','fname_mask');
    else
      load(fname_regmat);
    end;

    % resample fname_ref and fname_mask to atlas
    if ~parms.prereg_flag && parms.atlasflag
      % rigid body transform average to atlas
      %  and resample to 256^3, 1mm^3, coronal slices
      fprintf('%s: resampling %s to atlas and LIA...\n',mfilename,fname_ref);
      [vol,mr_parms] = ctx_load_mgh(fname_ref);
      vol = ctx_resample_to_LIA(vol,M_atl_to_subj);
      ctx_save_mgh(vol,fname_ref_res,mr_parms);
      % repeat for mask
      fprintf('%s: resampling %s to atlas and LIA...\n',mfilename,fname_mask);
      [volmask,mr_parms] = ctx_load_mgh(fname_mask);
      volmask = ctx_resample_to_LIA(volmask,M_atl_to_subj);
      volmask.imgs = 1.0*(volmask.imgs>0.5);
      ctx_save_mgh(volmask,fname_mask_res,mr_parms);
      clear vol volmask;
    elseif ~parms.atlasflag % create fname_ref_res in original space
      % copy fname_ref to fname_ref_res
      % it will overwrite fname_ref_res from registration step
      cmd = sprintf('cp %s %s',fname_ref,fname_ref_res);
      disp(cmd);
      [status,result] = unix(cmd);
      disp(result);
      if status, error('cmd %s failed',cmd); end;
      % copy fname_mask to fname_mask_res
      % it will overwrite fname_mask_res from registration step
      cmd = sprintf('cp %s %s',fname_mask,fname_mask_res);
      disp(cmd);
      [status,result] = unix(cmd);
      disp(result);
      if status, error('cmd %s failed',cmd); end;
      % set M_atl_to_subj 
      M_atl_to_subj = eye(4);
    end;
  end;

  % smooth brain mask for new reg
  if ~parms.oldreg_flag
    fname_mask_res_sm = sprintf('%s/%s_res_mask_sm%s',ContainerPath,outstem_mask,parms.ext); 
    if ~exist(fname_mask_res_sm,'file') || parms.forceflag
      mmil_dilate_mask([],'fname_in',fname_mask_res,'fname_out',fname_mask_res_sm,...
        'smooth1',parms.mask_smooth1,'thresh1',parms.mask_thresh1,...
        'smooth2',parms.mask_smooth2,'thresh2',parms.mask_thresh2,...
        'smooth3',parms.mask_smooth3,'forceflag',parms.forceflag);
    end;
  end;
else
  M_atl_to_subj = eye(4);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for t=1:length(parms.scantypes)
  fnamestem = parms.scantypes{t};
  switch fnamestem
    case {'FLASHhi','MEDIChi'}
      outstem = 'hiFA';
    case {'FLASHlo','MEDIClo'}
      outstem = 'loFA';
    case {'XetaT2'}
      outstem = 'T2w';
    otherwise
      outstem = fnamestem;
  end;

  if parms.wmbc_flag && ismember(fnamestem,parms.wmbc_scantypes)
    if ismember(fnamestem,parms.T2w_scantypes)
      image_range = parms.T2w_image_range;
    else
      image_range = parms.image_range;
    end;
  else
    image_range = [];
  end;
  
  % determine correct file name suffix
  suffix = set_suffix(parms,ScanInfo,fnamestem);
  if strcmp(suffix,'none'), continue; end;
  
  % create list of valid input files
  scanlist = {};
  tmpinfo = mmil_getfield(ScanInfo,fnamestem);
  for i = 1:length(tmpinfo)
    fname = sprintf('%s/%s%d%s%s',ContainerPath,fnamestem,i,suffix,parms.ext);
    if ~exist(fname,'file'), continue; end;
    if parms.rawQCflag &&...
       ((MPRflag && strcmp(fnamestem,'MPR')) ||...
        (FLASHhiflag && strcmp(fnamestem,'FLASHhi')))
      n = regexp({qcscans.fname},sprintf('%s%d',fnamestem,i)); n = [n{:}];
      if isempty(n)
        fprintf('%s: BAD QC WARNING: skipping file %s for bad QC\n',...
          mfilename,fname);
        fname = [];
      end;
    end;
    if ~isempty(fname), scanlist{end+1} = fname; end;
  end;

  % register to reference and average scans within type
  fname_avg = sprintf('%s/%s_res%s',...
    ContainerPath,outstem,parms.ext);     
  if strcmp(ref_scantype,fnamestem) && length(scanlist)==1
    if ~exist(fname_avg,'file') || parms.forceflag
      if parms.nativeflag || ~parms.atlasflag
        fname_ref = scanlist{1};
      elseif parms.atlasflag
        fname_ref = fname_ref_res;
      end;
      [vol,mr_parms] = ctx_load_mgh(fname_ref);
      if ~parms.nativeflag && ~parms.atlasflag
        % resample to 256^3, 1mm^3, coronal slices
        fprintf('%s: resampling %s to LIA...\n',mfilename,fname_ref);
        vol = ctx_resample_to_LIA(vol,eye(4));
      end;
      % clip values
      if ~isempty(image_range)
        vol.imgs(vol.imgs<image_range(1)) = image_range(1);
        vol.imgs(vol.imgs>image_range(2)) = image_range(2);
      end;
      ctx_save_mgh(vol,fname_avg,mr_parms);
    end;
  elseif length(scanlist)>=1
    if parms.oldreg_flag
      errcode = mmil_average_volumes_old(scanlist,fname_mask_res,...
        'fname_ref',fname_ref_res,...
        'fname_out',fname_avg,...
        'M_init',M_atl_to_subj,...
        'ref_scantype',ref_scantype,...
        'input_scantype',fnamestem,...
        'image_range',image_range,...
        'cleanupflag',parms.cleanupflag,...
        'forceflag',parms.forceflag);
    else
      errcode = mmil_average_volumes(scanlist,fname_mask_res_sm,...
        'fname_ref',fname_ref_res,...
        'fname_out',fname_avg,...
        'M_init',M_atl_to_subj,...
        'ref_scantype',ref_scantype,...
        'input_scantype',fnamestem,...
        'image_range',image_range,...
        'smoothmask_flag',false,...
        'cleanupflag',parms.cleanupflag,...
        'forceflag',parms.forceflag);
    end;  
  end;
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M_atl_to_subj_init = initial_rigidreg(fname,fname_out,M_offset,parms)
  
  [pathstr,fname_prefix,ext] = fileparts(fname_out);
  [pathstr_in,fname_in,ext_in] = fileparts(fname);

  % set output folder
  tmp_outdir = sprintf('%s/tmp_T1atl',pathstr);
  mmil_mkdir(tmp_outdir);

  % iteration 1
  tmp_outdir_itr1 = sprintf('%s/rr',tmp_outdir);
  fname_itr1 = sprintf('%s/%s_reg2atl_itr1.mgz',tmp_outdir_itr1,fname_in);
  M_atl_to_subj_itr1 = mmil_reg(parms.fname_atl,fname,'fname_maskA',parms.fname_atl_mask, ...
    'rigid_flag',1,'outdir',tmp_outdir_itr1);
  [vol_in,M_in] = fs_load_mgh(fname);
  [M_atl,volsz_atl] = fs_read_header(parms.fname_atl);
  [vol_itr1,M_itr1] = mmil_resample_vol(vol_in,M_in,...
    'M_ref',M_atl,'nvox_ref',volsz_atl(1:3),...
    'M_reg',M_atl_to_subj_itr1);
  clear vol_in;
  fs_save_mgh(vol_itr1,fname_itr1,M_itr1);
  clear vol_itr1;

  % iteration 2
  tmp_outdir_itr2 = sprintf('%s/rrjpdf',tmp_outdir);
  fname_itr2 = sprintf('%s/%s_reg2atl_itr2.mgz',tmp_outdir_itr2,fname_in);
  M_atl_to_subj_itr2 = mmil_reg(parms.fname_atl,fname_itr1,'fname_maskA',parms.fname_atl_maskBroad,...
    'rigid_flag',1,'options','-jpr -jpbrrMultScale','outdir',tmp_outdir_itr2);

  % remove temporary directory
  cmd = sprintf('rm -r %s',tmp_outdir);
  [status,result] = unix(cmd);
  if status
    error('failed to delete the temporary directory %s:\n%s\n',...
      tmp_outdir,result);
  end;
  
  M_atl_to_subj_init = (M_atl_to_subj_itr1 * M_atl_to_subj_itr2) * M_offset;
  
  % resampling to atlas with the offset
  M_LIA_to_orig = mmil_resample_to_LIA(fname,fname_out,'M_reg',M_atl_to_subj_init);
   
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function suffix = set_suffix(parms,ScanInfo,fnamestem)
  suffix = 'none';
  tmpinfo = mmil_getfield(ScanInfo,fnamestem);
  for i=1:length(tmpinfo)
    % check if grad unwarp should have been done
    unwarpflag = -1;
    if parms.gradunwarp_flag
      gradwarpinfo = mmil_getfield(tmpinfo(i),'gradwarpinfo');
      if ~isempty(gradwarpinfo)
        unwarpflag = mmil_getfield(gradwarpinfo,'unwarpflag',-1);
      end;
    end;
    if unwarpflag >= 0
      if parms.wmbc_flag && ismember(fnamestem,parms.wmbc_scantypes)
        suffixlist = {'_B1_uw_wmbc','_uw_wmbc'};
        if parms.gradunwarp_flag==1
          suffixlist = cat(2,suffixlist,{'_B1_wmbc','_wmbc'});
        end;
      elseif parms.nu_flag && ismember(fnamestem,parms.nu_scantypes)
        suffixlist = {'_B1_uw_nu','_uw_nu'};
        if parms.gradunwarp_flag==1
          suffixlist = cat(2,suffixlist,{'_B1_nu','_nu'});
        end;
      else
        suffixlist = {'_B1_uw','_uw'};
        if parms.gradunwarp_flag==1
          suffixlist = cat(2,suffixlist,{'_B1',''});
        end;
      end;
    else
      if parms.wmbc_flag && ismember(fnamestem,parms.wmbc_scantypes)
        suffixlist = {'_B1_wmbc','_wmbc'};
      elseif parms.nu_flag && ismember(fnamestem,parms.nu_scantypes)
        suffixlist = {'_B1_nu','_nu'};
      else
        suffixlist = {'_B1',''};
      end;
    end;
    for s=1:length(suffixlist)
      fname = sprintf('%s/%s%d%s%s',...
        parms.ContainerPath,fnamestem,i,suffixlist{s},parms.ext);
      if exist(fname,'file')
        suffix = suffixlist{s};
        break;
      end;
    end;
    if ~strcmp(suffix,'none'), break; end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

