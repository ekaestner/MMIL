function fs_mri_convert(fname_in,fname_out,varargin)
%function fs_mri_convert(fname_in,fname_out,[options])
%
% Purpose: call freesurfer's mri_convert for volume file format conversions
%
% Usage:
%  fs_mri_convert(fname_in,fname_out,'key1', value1,...); 
%
% Required Parameters:
%   fname_in: full or relative path name of input file
%   fname_out: full or relative path name of output file
%
% Optional parameters:
%  'out_orient': output orientation (e.g. 'LPS', 'RAS', etc.)
%     if empty or omitted, will keep original orientation
%    {default = []}
%  'options': other options (see mri_convert --help)
%    {default: []}
%  'verbose': [0|1] display status messages
%    {default: 0}
%  'forceflag': [0|1] overwrite existing output files
%    {default: 0}
%  'TR': repetition time  
%    {default: []}
% Created:  07/19/08 Don Hagler
% Prev Mod: 06/16/14 Don Hagler
% Last Mod: 08/21/13 Dani Cornejo
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin, 2), return; end;
parms = mmil_args2parms(varargin, { ...
  'out_orient',[],[],...
  'options',[],[],...
  'verbose',false,[false true],...
  'forceflag',false,[false true],...
  'TR',[],[0 Inf],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(fname_in,fname_out)
  error('fname_in and fname_out must be different');
end;
if ~exist(fname_in,'file')
  error('file %s not found',fname_in);
end;

if ~exist(fname_out,'file') || parms.forceflag
  if parms.verbose
    fprintf('%s: converting %s to %s...\n',mfilename,fname_in,fname_out);
  end;
  cmd = 'mri_convert';
  if ~isempty(parms.out_orient)
    cmd = [cmd ' --out_orientation ' upper(parms.out_orient)];
  end;
  if ~isempty(parms.TR)
    % checking TR units 
    if floor(parms.TR) ~= parms.TR
      TR = num2str(floor(parms.TR*1000));     
      cmd = [cmd ' -tr ' TR];
    else 
      TR = num2str(floor(parms.TR));     
      cmd = [cmd ' -tr ' TR];  
    end 
  end;
  if ~isempty(parms.options)
    cmd = [cmd ' ' parms.options];
  end;
  cmd = [cmd ' ' fname_in ' ' fname_out]; 
  [status,msg] = unix(cmd);
  if status || ~exist(fname_out,'file')
    error('failed to convert %s to %s:\n%s\n%s\n',...
      fname_in,fname_out,cmd,msg);
  end;
end;

