function [ContainerOutPath,errcode] = MMIL_Import_MEG(ContainerPath,RootDirs,varargin)
%function [ContainerOutPath,errcode] = MMIL_Import_MEG(ContainerPath,RootDirs,varargin)
%
% Required Input:
%  ContainerPath: full path name of directory containing original MEG data
%  RootDirs:
%    a struct which must contain the following fields:
%         orig_meg, raw_meg, proc_meg
%    these specify the locations of data
%   orig_meg: full path root dir for orig MEG Containers
%   raw_meg: full path root dir for raw MEG Containers
%   proc_meg: full path root dir for processed MEG Containers
%
% Optional Input:
%  'rawflag': [0|1|2] copy data from orig to raw (no maxfilter)
%    0: run maxfilter, even if no corrections performed
%       (may be necessary for some versions of data)
%    1: do not run maxfilter, just link data from orig to raw
%    2: do not run maxfilter, just copy data from orig to raw
%    ignored if sssflag = 1 or moveflag = 1
%    {default = 0}
%  'sssflag': [0|1] run maxfilter with signal space separation
%    {default = 0}
%  'moveflag': [0|1] run maxfilter with maxmove
%    Requires continuous HPI coil measurements
%    If 1, ssflag is ignored
%    {default = 0}
%  'format': output data format (short, long, or float)
%     {default = 'float'}
%  'forceflag': [0|1] whether to overwrite existing output
%    {default = 0}
%
% Output:
%   ContainerOutDir: name of output MEGRAW Container
%     (after running maxfilter)
%   errcode: 0 if no errors, 1 if errors
%
% Created:  02/17/11 by Don Hagler
% Last Mod: 08/08/16 by Don Hagler
%

ContainerOutDir = [];
errcode = 0;

if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin,{...
  'rawflag',0,[0:2],...
  'sssflag',false,[false true],...
  'moveflag',false,[false true],...
  'format','float',{'short','long','float'},...
  'forceflag',false,[false true],...
...
  'st',6,[0.1 100],...
  'corr',0.95,[0 1],...
...
  'exclude_files',{'avg','onlineavg','emptyroom','resting'},[],...
});

fprintf('%s: Importing %s...\n',mfilename,ContainerPath);

[tmp,ContainerDir,text] = fileparts(ContainerPath);
if ~isempty(text)
  ContainerDir = [ContainerDir text];
end;

% replace any '.' in VisitID with '_'
%  (causes problems when inserted in job names and for ADNI)
VisitID = regexprep(ContainerDir,'\.','_');

% input file list
flist = dir(sprintf('%s/*.fif',ContainerPath));
nfiles = length(flist);
if ~nfiles
  fprintf('%s: ERROR: no fif files found in %s\n',...
    mfilename,ContainerPath);
  errcode = 1;
  return;
end;
filenames = {flist.name};

% exclude useless files that will cause problems for processing
ind_exclude = [];
for f=1:nfiles
  for x=1:length(parms.exclude_files)
    if ~isempty(regexp(filenames{f},parms.exclude_files{x}))
      ind_exclude = [ind_exclude,f];
    end;
  end;
end;
ind_keep = setdiff(1:nfiles,ind_exclude);
filenames = filenames(ind_keep);
nfiles = length(filenames);
if ~nfiles
  fprintf('%s: ERROR: no valid fif files found in %s\n',...
    mfilename,ContainerPath);
  errcode = 1;
  return;
end;

% create output directory
ContainerOutPath = [RootDirs.raw_meg '/MEGRAW_' VisitID];
[succ,msg] = mkdir(ContainerOutPath);

% run maxfilter for each file
input_data_files = cell(nfiles,1);
output_data_files = cell(nfiles,1);
fname_ref = [ContainerPath '/' filenames{1}];
k = 1;
for f=1:nfiles
  fname_in = [ContainerPath '/' filenames{f}];
  fname_out = [ContainerOutPath '/' filenames{f}];
  if parms.moveflag
    fname_out = regexprep(fname_out,'.fif','-stmc.fif');
    fname_log = regexprep(fname_out,'.fif','-stmc.log');
  elseif parms.sssflag
    fname_out = regexprep(fname_out,'.fif','-sss.fif');
  end;
  if ~exist(fname_out,'file') || parms.forceflag
    if parms.rawflag && ~parms.sssflag && ~parms.moveflag
      if parms.rawflag == 1
        cmd = sprintf('ln -s %s %s',fname_in,fname_out);
      elseif parms.rawflag == 2
        cmd = sprintf('cp -p %s %s',fname_in,fname_out);
      end;
      fprintf('%s: cmd = \n%s\n\n',mfilename,cmd);
      [s,r] = mmil_unix(cmd);
    else
      cmd = 'nice maxfilter';
      cmd = sprintf('%s -f %s -o %s',cmd,fname_in,fname_out);
      cmd = sprintf('%s -v -format %s',cmd,parms.format);
      if parms.sssflag || parms.moveflag
        cmd = sprintf('%s -st %0.2f',cmd,parms.st);
        cmd = sprintf('%s -corr %0.2f',cmd,parms.corr);
      end;
      if parms.moveflag
        cmd = sprintf('%s -movecomp',cmd);
        if f>1
          cmd = sprintf('%s -trans %s',cmd,fname_ref);
        end;      
        cmd = sprintf('%s -hpicons -force -autobad on | tee -a %s',cmd,fname_log);
      end;
      fprintf('%s: cmd = \n%s\n\n',mfilename,cmd);
      [s,r] = mmil_unix(cmd);
      if s && ~isempty(regexp(r,'SSS was already applied'))
        fprintf('%s: WARNING:\n%s\n',mfilename,r);
        cmd = sprintf('ln -s %s %s',fname_in,fname_out);
        fprintf('%s: cmd = \n%s\n\n',mfilename,cmd);
        [s,r] = mmil_unix(cmd);
      end;
    end;
    if s
      fprintf('%s: WARNING: failed to import raw MEG data:\n%s\n',...
        mfilename,r);
      errcode = 1;
    else
      fprintf('%s\n\n',r);
    end;
  end;
  if exist(fname_out,'file')
    input_data_files{k} = fname_in;
    output_data_files{k} = fname_out;
    k = k + 1;
  end;
end;

% create ContainerInfo.mat
ContainerInfo = [];
ContainerInfo.ContainerType = 'MEGRAW';
ContainerInfo.VisitID = VisitID;
ContainerInfo.input_data_files = input_data_files;
ContainerInfo.output_data_files = output_data_files;
errcode = MMIL_Save_ContainerInfo(ContainerOutPath,ContainerInfo);



