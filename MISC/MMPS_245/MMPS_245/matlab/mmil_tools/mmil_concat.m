function mmil_concat(fnamelist,fname,varargin)
%function mmil_concat(fnamelist,fname,[options])
%
% Purpose: concatenate a series of mgh/mgz files
%   this is a wrapper for FreeSurfer's mri_concat
%
% Required Input:
%   fnamelist: cell array of full path file names
%   fname: output file name
%
% Optional Parameters:
%   'meanflag': [0|1] calculate mean instead of concatenating
%     {default = 0}
%   'options': option string to use any of the mri_concat command line options
%     {default = []}
%   'forceflag': overwrite existing output
%     {default = 0}
%
% Created:  10/30/12 by Don Hagler
% Prev Mod: 08/06/14 by Don Hagler
% Last Mod: 11/19/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin,{,...
  'meanflag',false,[false true],...
  'options',[],[],...
  'forceflag',false,[false true],...
...
  'check_size_flag',true,[false true],...
  'cmd_name','mri_concat_530',[],...
});

if exist(fname,'file') && ~parms.forceflag, return; end;

if ~iscell(fnamelist), fnamelist = {fnamelist}; end;
nfiles = length(fnamelist);
for f=1:nfiles
  if ~exist(fnamelist{f},'file')
    error('file %s not found',fnamelist{f});
  end;
end;

% get output directory
if mmil_isrelative(fname), fname = [pwd '/' fname]; end;
[tpath,tstem] = fileparts(fname);
mmil_mkdir(tpath);
fname_list = [tpath '/' tstem '_fname_list.txt'];

% set temporary output name for mean (to allow reshaping output)
if parms.check_size_flag
  [tmp,tmp_infix] = fileparts(tempname);
  fname_tmp = [tpath '/' tstem '_' tmp_infix '.mgh'];
else
  fname_tmp = fname;
end;

% write temporary file with list of input files
fid = fopen(fname_list,'w');
if fid==-1
  error('failed to open %s for writing',fname_list);
end;
fprintf(fid,'%s\n',fnamelist{:});
fclose(fid);

% make unix call
cmd = sprintf('%s --f %s --o %s',parms.cmd_name,fname_list,fname_tmp);
if parms.meanflag
  cmd = [cmd ' --mean'];
end;
if ~isempty(parms.options)
  cmd = [cmd ' ' parms.options];
end;
[s,r] = mmil_unix(cmd);
if s, error('command "%s" failed:\n%s',cmd,r); end;

if parms.check_size_flag
  [vals,M,tmp,volsz] = fs_load_mgh(fname_tmp);
  % check for surface file with wrong order of dimensions
  if (length(find(volsz(1:3)==1)) == 2) && volsz(1)==1
    vals = reshape(vals,[max(volsz(1:3)),1,1,volsz(4)]);
  end;
  fs_save_mgh(vals,fname,M);
  % remove temporary file
  delete(fname_tmp);
end;

% remove file list file
delete(fname_list);

