function status=fs_warp_vol2atlas(subj,fname_in,fname_out,varargin);
%function status=fs_warp_vol2atlas(subj,fname_in,fname_out,[options]);
%
% Usage:
%  [status] = fs_warp_vol2atlas(subj,fname_in,fname_out,'key1', value1,...);
%
% Required input:
%  subj: string specifying the subject name
%  fname_in: full or relative path of input file name
%  fname_out: full or relative path of output file name
%
% Optional parameters:
%  'transform' - non-linear transform filename
%    {default: 'subjdir/subj/mri/transforms/talairach.m3z'}
%  'inverse_flag' - [0|1] toggle apply inverse transform (warp atlas to subject)
%    {default: 0}
%  'resamp_type' - resample type
%    ('interpolate', 'weighted', 'nearest', 'sinc', 'cubic')
%    {default: 'nearest'}
%  'subjdir' - subjects directory (override SUBJECTS_DIR environment variable)
%    subjdir/subj should contain the freesurfer subject directory
%    {default = $SUBJECTS_DIR}
%  'overwrite_flag' - [0|1] toggle overwrite existing output files
%    {default: 0}
%
% created: 11/17/06 by Don Hagler
% last modified: 02/08/07 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse options

status = 0;

try
  options = varargin;
  for index = 1:length(options)
      if iscell(options{index}) & ~iscell(options{index}{1}), options{index} = { options{index} }; end;
  end;
  if ~isempty( varargin ), g=struct(options{:}); 
  else g = []; end;
catch
  fprintf('%s: ERROR: bad input arguments: use {''key'', value, ... }\n',mfilename);
  return;
end;

if ~isfield(g,'transform'), g.transform = []; end;
if ~isfield(g,'resamp_type'), g.resamp_type = 'nearest'; end;
if ~isfield(g,'subjdir'), g.subjdir = []; end;
if ~isfield(g,'inverse_flag'), g.inverse_flag = 0; end;
if ~isfield(g,'overwrite_flag'), g.overwrite_flag = 0; end;

gfields = fieldnames(g);
for index=1:length(gfields)
   switch gfields{index}
   case {'transform' 'resamp_type'...
         'subjdir' 'overwrite_flag' 'inverse_flag'},;
   otherwise, error([mfilename ': unrecognized option: ''' gfields{index} '''' ]);
   end;
end;

% get rid of options struct
transform = g.transform;
resamp_type = g.resamp_type;
subjdir = g.subjdir;
inverse_flag = g.inverse_flag;
overwrite_flag = g.overwrite_flag;
clear g;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check parameters
if nargin<3, help(mfilename); return; end;  

if ~exist('subjdir','var'), subjdir = []; end;
if isempty(subjdir)
  subjdir = getenv('SUBJECTS_DIR');
  if isempty(subjdir)
    fprintf('%s: SUBJECTS_DIR not defined as an environment variable... quitting\n',mfilename);
    return;
  end;
end;
setenv('SUBJECTS_DIR',subjdir);

if ~exist(fname_in,'file')
  fprintf('%s: ERROR: input file %s not found\n',mfilename,fname_in);
  return;
end;

if isempty(transform)
  transform = sprintf('%s/%s/mri/transforms/talairach.m3z',subjdir,subj);
end;
if ~exist(transform,'file')
  fprintf('%s: transform file %s not found...quitting\n',mfilename,transform);
  return;
end;

if ~ismember(resamp_type,{'interpolate' 'weighted' 'nearest' 'sinc' 'cubic'})
  fprintf('%s: unsupported resamp_type: %s...quitting\n',mfilename,resamp_type);
  return;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(fname_out,'file') | overwrite_flag
  % create mri_convert cmd string
  if ~inverse_flag
    cmd = sprintf('mri_convert -rt %s -at %s %s %s',...
      resamp_type,transform,fname_in,fname_out);
  else
    cmd = sprintf('mri_convert -rt %s -ait %s %s %s',...
      resamp_type,transform,fname_in,fname_out);
  end;
  % run cmd
  fprintf('%s\n',cmd);
  [status,result]=unix(cmd);
  if status
    fprintf('%s: ERROR:\n',mfilename);
    disp(result);
  end;
end;
