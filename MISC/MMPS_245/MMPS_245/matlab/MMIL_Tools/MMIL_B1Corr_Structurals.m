function errcode = MMIL_B1Corr_Structurals(ContainerPath,varargin)
%function errcode = MMIL_B1Corr_Structurals(ContainerPath,[options])
%
% Purpose: correct structural images for B1 intensity inhomogeneities
%
% Required Input:
%   ContainerPath: full path of MRIPROC container
%
% Optional Input:
%  'scantypes': cell array of scan types to correct
%    {default = {'MPR','FLASHhi','FLASHlo','MEDIChi','MEDIClo','XetaT2'}}
%  'ext': file extension (i.e. '.mgh' or '.mgz')
%    {default = '.mgz'}
%  'forceflag': [0|1] whether to overwrite existing output
%    {default = 0}
%
% Created:  04/12/12 by Vijay Venkatraman
% Last Mod: 04/17/14 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errcode = 0;
if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'scantypes',{'MPR','FLASHhi','FLASHlo','MEDIChi','MEDIClo','XetaT2'},[],...
  'ext','.mgz',{'.mgh','.mgz'},...  
  'forceflag',false,[false true],...
... % hidden
  'scantypes_GEcal',{'MPR','FLASHhi','FLASHlo'},[],...
  'manufacturers',{'GE MEDICAL SYSTEMS','SIEMENS','PHILIPS MEDICAL SYSTEMS'},[],...
  'manufacturers_corr',{'GE MEDICAL SYSTEMS','SIEMENS'},[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load container info
[ContainerInfo,errcode] = MMIL_Load_ContainerInfo(ContainerPath);
if errcode ~= 0, return; end;

if isempty(ContainerInfo)
  fprintf('ERROR: ContainerInfo is empty\n');
  errcode = 1;
  return;
end

if ~isfield(ContainerInfo,'ScanInfo')
  fprintf('ERROR: ContainerInfo.ScanInfo missing\n');
  errcode = 1;
  return;
end
ScanInfo = ContainerInfo.ScanInfo;

% check manufacturer
parms.manufacturer = strtrim(upper(ContainerInfo.Manufacturer));
if ~ismember(parms.manufacturer,parms.manufacturers)
  fprintf('%s: ERROR: Manufacturer %s not supported\n',...
    mfilename,parms.manufacturer);
  errcode = 1;
  return;
end;
if ~ismember(parms.manufacturer,parms.manufacturers_corr)
  fprintf('%s: B1 correction not required for Manufacturer %s...\n',...
    mfilename,parms.manufacturer);
  return;
end;

% check if B1 correction is needed (for at least one scan)
parms.B1corr_flag = mmil_check_B1corr(ScanInfo,parms.scantypes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate B1 correction field
ctx_volrat = [];
if parms.B1corr_flag
  % calculate B1 correction field from GE calibration scan
  fname = sprintf('%s/GEB1CAL_ratio%s',ContainerPath,parms.ext);
  if ~exist(fname,'file') || parms.forceflag
    % select high resolution scan to register B1cal scan to  
    stype = [];
    fstem_cal = 'GEB1CAL1';
    if exist(sprintf('%s/%s%s',ContainerPath,fstem_cal,parms.ext),'file')
      for f=1:length(parms.scantypes_GEcal)
        fnamestem = parms.scantypes_GEcal{f};
        fstem_hires = sprintf('%s1',fnamestem);
        fname_hires = sprintf('%s/%s%s',ContainerPath,fstem_hires,parms.ext);
        if exist(fname_hires,'file')
          % check whether B1 correction needed for this scantype
          if mmil_check_B1corr(ScanInfo,fnamestem)
            stype = fnamestem;
            break;
          end;
        end;
        if ~isempty(stype), break; end;
      end;
      if ~isempty(stype)
        [ctx_volrat,mr_parms] = ...
          mmil_estimate_B1RxRatio_GECal(ContainerPath,fstem_cal,...
          'fstem_hires',fstem_hires,'stype',stype,'ext',parms.ext);
      end;
      if ~isempty(ctx_volrat)
        ctx_save_mgh(ctx_volrat,fname,mr_parms);
      end;
    end;
  elseif exist(fname,'file')
    [ctx_volrat,mr_parms] = ctx_load_mgh(fname);
  end;

  % calculate B1 correction field from ADNI calibration scans
  if isempty(ctx_volrat)
    fname = sprintf('%s/B1volrat%s',ContainerPath,parms.ext);
    if ~exist(fname,'file') || parms.forceflag
      BCmstem = 'B1BC_mag1';
      BCpstem = 'B1BC_phi2';
      HCmstem = 'B1HC_mag1';
      HCpstem = 'B1HC_phi2';
      ctx_volrat = mmil_estimate_B1RxRatio_ADNI(ContainerPath,...
        BCmstem,BCpstem,HCmstem,HCpstem,parms.ext);
      if ~isempty(ctx_volrat)
        ctx_save_mgh(ctx_volrat,fname);
      end;
    else
      [ctx_volrat,mr_parms] = ctx_load_mgh(fname);
    end;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% apply corrections for each scan type
for f=1:length(parms.scantypes)
  fnamestem = parms.scantypes{f};
  tmpinfo_array = mmil_getfield(ScanInfo,fnamestem,[]);
  for i = 1:length(tmpinfo_array)
    % B1 corrections
    if parms.B1corr_flag
      fname_in = sprintf('%s/%s%d%s',ContainerPath,fnamestem,i,parms.ext);
      fname_out = sprintf('%s/%s%d_B1%s',ContainerPath,fnamestem,i,parms.ext);
      % check whether B1 correction needed for this scan
      if mmil_check_B1corr(ScanInfo,fnamestem,i) && (~isempty(ctx_volrat) && ...
          exist(fname_in,'file')) && (~exist(fname_out,'file') || parms.forceflag)
        fprintf('%s: applying B1 correction for %s...\n',mfilename,fname_in); 
        [ctx_vol,mr_parms] = ctx_load_mgh(fname_in);
        ctx_volrat_res = vol_resample(ctx_volrat, ctx_vol, eye(4), 1);
        ctx_vol.imgs = ctx_vol.imgs.*ctx_volrat_res.imgs;
        ctx_save_mgh(ctx_vol,fname_out,mr_parms);
      end;
    end;
  end;
end;

return;
