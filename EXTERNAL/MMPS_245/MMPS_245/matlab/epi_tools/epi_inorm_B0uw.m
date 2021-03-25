function fname_rev_inorm = epi_inorm_B0uw(fname_for,fname_rev,varargin)
%function fname_rev_inorm = epi_inorm_B0uw(fname_for,fname_rev,varargin)
%
% Required Parameters:
%   fname_for: name of input image with "forward" phase-encode polarity
%   fname_rev: name of input image with "reverse" phase-encode polarity
%
% Optional Parameters:
%  'fname_mask': name of mask image used for global intensity normalization.
%     if empty, will generate mask from forward image
%     {default = []} 
%  'fname_input_param': filename of parameter file used by reg binary
%     {default = []}
%  'fname_rev_inorm': output filename of intensity normalized "reverse"
%     phase-encode polarity; if empty, generated from fname_rev
%     {default = []}
%  'cleanup_flag': [0|1] delete mask file images created if needed 
%     {default = 0}
%  'forceflag': [0|1] whether to overwrite existing output files
%     {default = 0}
%
% Created : 07/27/12 by Vijay Venkatraman
% Last Mod: 10/27/12 by Don Hagler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;

parms = mmil_args2parms(varargin,{...
  'fname_mask',[],[],...
  'fname_rev_inorm',[],[],...
  'fname_input_param',[],[],...
  'cleanup_flag',false,[false true],...
  'forceflag',false,[false true],...
});

if ~exist(fname_for,'file')
  error('input file %s not found',fname_for);
end;
if ~exist(fname_rev,'file')
  error('input file %s not found',fname_rev);
end;
if isempty(parms.fname_mask)
  [for_pathstr,for_name,for_ext] = fileparts(fname_for);
  parms.fname_mask = sprintf('%s/%s_mask%s',for_pathstr,for_name,for_ext);
  if ~exist(parms.fname_mask,'file')
    epi_brainmask(fname_for,parms.fname_mask,'forceflag',parms.forceflag);
    parms.cleanup_flag = 1;
  end;
end;

% set output name
[pathstr, name, ext] = fileparts(fname_rev);
if isempty(parms.fname_rev_inorm)
  parms.fname_rev_inorm = sprintf('%s/%s_inorm%s',pathstr,name,ext);
end;

if parms.forceflag || ~exist(parms.fname_rev_inorm,'file')
  % write the input parameter file
  if isempty(parms.fname_input_param)
    parms.fname_input_param = sprintf('%s/inputParamsScale.txt',pathstr);
    fid=fopen(parms.fname_input_param,'wt');
    if fid==-1
      error('failed to open file %s for writing',parms.fname_input_param);
    end;
    fprintf(fid,'standardizeIntensities_S_To_T_Only = 1 \n');
    fprintf(fid,'rescaleToTarget = 1  Voxel dims: override defaults \n');
  end;
  % global intensity normalize using 'reg' command
  cmd = sprintf('reg -t %s -s %s -tm %s -ip %s -od %s',fname_for,...
    fname_rev,parms.fname_mask,parms.fname_input_param,pathstr);
  [status, result] = unix(cmd);
  if status
    error('%s: %s failed, reason: %s \n',mfilename,cmd,result);
  end;
  % rename the output file
  fname_rev_out = sprintf('%s/%s_GlobalIntensityNorm%s',pathstr,name,ext);
  fname_inorm_txt = sprintf('%s/globalIScale.txt',pathstr);
  cmd = sprintf('mv %s %s',fname_rev_out,parms.fname_rev_inorm);
  [status,result] = unix(cmd);
  if status
    error('%s: %s failed, reason: %s \n',mfilename,cmd,result);
  end;
  % clean up output text file from global intensity normalization
  cmd = sprintf('rm %s',fname_inorm_txt);
  [status,result] = unix(cmd);
  if status
    error('%s: %s failed, reason: %s \n',mfilename,cmd,result);
  end;
  % clean up input parameter file
  cmd = sprintf('rm %s',parms.fname_input_param);
  [status,result] = unix(cmd);
  if status
    error('%s: %s failed, reason: %s \n',mfilename,cmd,result);
  end;
end;

if parms.cleanup_flag
  cmd = sprintf('rm %s',parms.fname_mask);
  [status,result] = unix(cmd);
  if status
    error('%s: %s failed, reason: %s \n',mfilename,cmd,result);
  end;
end;

fname_rev_inorm = parms.fname_rev_inorm;

return;





