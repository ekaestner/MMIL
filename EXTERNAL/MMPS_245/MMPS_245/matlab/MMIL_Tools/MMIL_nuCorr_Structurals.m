function errcode = MMIL_nuCorr_Structurals(ContainerPath,varargin)
%function errcode = MMIL_nuCorr_Structurals(ContainerPath,[options])
%
% Purpose: correct structural images for bias-field intensity inhomogeneities
%
% Required Input:
%   ContainerPath: full path of MRIPROC container
%
% Optional Input:
%  'scantypes': cell array of scan types to correct
%     {default = {'MPR','FLASHhi','MEDIChi',}}
%  'ext': file extension (i.e. '.mgh' or '.mgz')
%    {default = '.mgz'}
%  'tal_flag': [0|1] register to Talairach space
%    for standardization of white matter intensity values
%    {default = 0}
%  'forceflag': [0|1] whether to overwrite existing output
%    {default = 0}
%
% Created:  04/12/12 by Vijay Venkatraman
% Last Mod: 11/27/15 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errcode = 0;
if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'scantypes',{'MPR','FLASHhi','MEDIChi'},[],...
  'ext','.mgz',{'.mgh','.mgz'},...
  'tal_flag',false,[false true],...
  'forceflag',false,[false true],...
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

% suffixlist allows to check whether input file exist
% checks for the input file in the order specified
suffixlist = {'_B1_uw','B1','_uw',''};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% apply nu corrections
for f=1:length(parms.scantypes)
  fnamestem = parms.scantypes{f};
  tmpinfo_array = mmil_getfield(ScanInfo,fnamestem,[]);
  for i = 1:length(tmpinfo_array)
    for s = 1:length(suffixlist)
      fname_in = sprintf('%s/%s%d%s%s',...
        ContainerPath,fnamestem,i,suffixlist{s},parms.ext);
      fname_out = sprintf('%s/%s%d%s_nu%s',...
        ContainerPath,fnamestem,i,suffixlist{s},parms.ext);
      if exist(fname_in, 'file')
        if ~exist(fname_out,'file') || parms.forceflag
          fprintf('%s: nu correction %s...\n',mfilename,fname_in);
          fs_nu_corr(fname_in,'fname_out',fname_out,...
            'tal_flag',parms.tal_flag,...
            'tmpdir',[ContainerPath '/tmp_nu_corr'],...
            'forceflag',parms.forceflag);
        end;
        break; % breaks from suffixlist loop after first fname_in occurrence
      end;
    end;  
  end;
end;

return;

