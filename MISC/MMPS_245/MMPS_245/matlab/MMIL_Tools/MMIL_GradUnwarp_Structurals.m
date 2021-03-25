function errcode = MMIL_GradUnwarp_Structurals(ContainerPath,varargin)
%function errcode = MMIL_GradUnwarp_Structurals(ContainerPath,[options])
%
% Purpose: correct structural images for grad warp distortions
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
% Last Mod: 11/15/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errcode = 0;
if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'scantypes',{'MPR','FLASHhi','FLASHlo','MEDIChi','MEDIClo','XetaT2'},[],...
  'ext','.mgz',{'.mgh','.mgz'},...  
  'forceflag',false,[false true],...
... % hidden
  'scantypes_medic',{'MEDIChi','MEDIClo'},[],...
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
  fprintf('%s: ERROR: manufacturer %s not supported\n',...
    mfilename,parms.manufacturer);
  errcode = 1;
  return;
end;
if ~ismember(parms.manufacturer,parms.manufacturers_corr)
  fprintf('%s: grad unwarp not required for manufacturer %s...\n',...
    mfilename,parms.manufacturer);
  return;
end;

% suffixlist allows to check whether B1 corrected file exist
% checks for the input file in the order specified
suffixlist = {'_B1',''};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% apply corrections for each scan type
for f=1:length(parms.scantypes)
  fnamestem = parms.scantypes{f};
  tmpinfo_array = mmil_getfield(ScanInfo,fnamestem,[]);
  for i=1:length(tmpinfo_array)
    tmpinfo = tmpinfo_array(i);
    gradwarpinfo = tmpinfo.gradwarpinfo;
    for s=1:length(suffixlist)
      fname_in = sprintf('%s/%s%d%s%s',...
        ContainerPath,fnamestem,i,suffixlist{s},parms.ext);
      if exist(fname_in,'file')
        fname_out = sprintf('%s/%s%d%s_uw%s',...
          ContainerPath,fnamestem,i,suffixlist{s},parms.ext);
        errcode = apply_gradunwarp(fname_in,fname_out,...
          gradwarpinfo,parms.forceflag);
        if errcode, return; end;
        break; % apply correction only to first in sufflix list
      end;  
    end;
    if length(tmpinfo.TE) > 1
      for j=1:length(tmpinfo.TE)
        fname_in = sprintf('%s/%s%d_%d%s',...
          ContainerPath,fnamestem,i,j,parms.ext);
        fname_out = sprintf('%s/%s%d_%d_uw%s',...
          ContainerPath,fnamestem,i,j,parms.ext);
        if exist(fname_in,'file')
          errcode = apply_gradunwarp(fname_in,fname_out,...
            gradwarpinfo,parms.forceflag);
          if errcode, return; end;
        end;
      end;
    end;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grad warp correction for medic R2 images
for f=1:length(parms.scantypes_medic)
  fnamestem = parms.scantypes_medic{f};
  tmpinfo_array = mmil_getfield(ScanInfo,fnamestem,[]);
  for i=1:length(tmpinfo_array)
    tmpinfo = tmpinfo_array(i);
    fname_in = sprintf('%s/%s%d_R2%s',ContainerPath,fnamestem,i,parms.ext);
    fname_out = sprintf('%s/%s%d_R2_uw%s',ContainerPath,fnamestem,i,parms.ext);
    if exist(fname_in,'file')
      gradwarpinfo = tmpinfo.gradwarpinfo;
      if isempty(gradwarpinfo) || ...
          ~isfield(gradwarpinfo,'gwtype') || isempty(gradwarpinfo.gwtype)
        fprintf('%s: ERROR: missing gradwarpinfo for %s\n',mfilename,fname_in);
        errcode = 1;
        return;
      elseif ~(exist(fname_out,'file')) || parms.forceflag
        fprintf('%s: unwarping %s...\n',mfilename,fname_in);
        [ctx_vol,mr_parms] = ctx_load_mgh(fname_in);
        ctx_voluw = ctx_unwarp_grad2(ctx_vol,...
          gradwarpinfo.gwtype,gradwarpinfo.unwarpflag,gradwarpinfo.isoctrflag,0);
        ctx_save_mgh(ctx_voluw,fname_out,mr_parms);
      end;
    end;
  end;
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errcode = apply_gradunwarp(fname_in,fname_out,gradwarpinfo,forceflag)
  errcode = 0;
  if isempty(gradwarpinfo) || ~isfield(gradwarpinfo,'gwtype') || ...
      isempty(gradwarpinfo.gwtype)
    fprintf('%s: ERROR: missing gradwarpinfo for %s\n',...
      mfilename,fname_in);
    errcode = 1;
  elseif ~exist(fname_out,'file') || forceflag
    fprintf('%s: unwarping %s...\n',mfilename,fname_in);
    [ctx_vol,mr_parms] = ctx_load_mgh(fname_in);
    ctx_voluw = ctx_unwarp_grad(ctx_vol,gradwarpinfo.gwtype,...
      gradwarpinfo.unwarpflag,gradwarpinfo.isoctrflag);
    ctx_save_mgh(ctx_voluw,fname_out,mr_parms);
  end;
return;
