function errcode = MMIL_Get_Gradwarp_Signature(ContainerPath,varargin)
%function errcode = MMIL_Get_Gradwarp_Signature(ContainerPath,[options])
%
% Purpose: determine gradwarp type for scans with ambiguous or missing header info
%
% Required Input:
%   ContainerPath: full path of MRIPROC container
%
% Optional Input:
%  'gradunwarp_flag' - [0|1] whether to correct for gradient warping
%     If 0, this function will return immediately
%     {default = 1}
%  'scantypes': cell array of scan types for which to find gradwarp signature
%     {default={'MPR','FLASHhi','FLASHlo','MEDIChi','MEDIClo','XetaT2'}}
%  'ext': file extension (i.e. '.mgh' or '.mgz')
%    {default = '.mgz'}
%  'forceflag': [0|1] whether to overwrite existing output
%    {default = 0}
%
% Created:  11/16/10 by Don Hagler
% Last Mod: 02/12/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errcode = 0;
if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'gradunwarp_flag',true,[false true],...  
  'scantypes',{'MPR','FLASHhi','FLASHlo','MEDIChi','MEDIClo','XetaT2'},[],...
  'ext','.mgz',{'.mgh','.mgz'},...
  'forceflag',false,[false true],...
... % hidden
  'manufacturers_corr',{'GE MEDICAL SYSTEMS'},[],... % may have ambiguous type
  'gwtypelist',[2,3,8],[],... % possible GE gw types
  'min_corr',0.9,[0.1,0.99],... % minimum correlation to be unambiguous
  'thresh',1,[0,1e6],... % threshold applied to image to detect edge
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~parms.gradunwarp_flag, return; end;

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
if ~ismember(parms.manufacturer,parms.manufacturers_corr)
  return;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_out = sprintf('%s/ScanInfo_GWS.mat',ContainerPath);
if ~exist(fname_out,'file')
  save_flag = 0;
  for f=1:length(parms.scantypes)
    fnamestem = parms.scantypes{f};
    tmpinfo = mmil_getfield(ScanInfo,fnamestem,[]);
    for i = 1:length(tmpinfo)
      fname = sprintf('%s/%s%d%s',ContainerPath,fnamestem,i,parms.ext);
      vol = [];
      gradwarpinfo = mmil_getfield(tmpinfo(i),'gradwarpinfo',[]);
      if isempty(gradwarpinfo), continue; end;
      unwarpflag = mmil_getfield(gradwarpinfo,'unwarpflag',0);
      if unwarpflag && isfield(gradwarpinfo,'ambiguousgwtype') &&...
        (~isfield(gradwarpinfo,'gwtype') || parms.forceflag)
        fprintf('%s: finding gradwarp type for %s...\n',mfilename,fname);
        vol = ctx_load_mgh(fname);
        if ndims(vol.imgs)<3
          fprintf('%s: WARNING: image is not a volume!\n',mfilename);
          return;
        end;
        [gradwarpinfo.gwtype,gradwarpinfo.gwcorr,gradwarpinfo.gwtypelist] =...
          ctx_get_gwtype(vol,...
            'isoctrflag',gradwarpinfo.isoctrflag,'gwtypelist',parms.gwtypelist,...
            'min_corr',parms.min_corr,'thresh',parms.thresh);
        if ~isempty(gradwarpinfo.gwtype)
          gradwarpinfo.ambiguousgwtype = 0;
        end;
        ScanInfo.(fnamestem)(i).gradwarpinfo = gradwarpinfo;
        save_flag = 1;
      end;
    end;
  end;
  if save_flag
    save(fname_out,'ScanInfo');
  end;
end;

if exist(fname_out,'file')
  load(fname_out);
  ContainerInfo.ScanInfo = ScanInfo;
  errcode = MMIL_Save_ContainerInfo(ContainerPath,ContainerInfo);
end;

