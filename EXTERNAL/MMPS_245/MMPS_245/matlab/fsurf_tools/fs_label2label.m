function fname_out = fs_label2label(fname,varargin)
%function fname_out = fs_label2label(fname,[options])
%
% Purpose: resample a label from one subject to another
%   (wrapper for freesurfer's mri_label2label)
%
% Usage:
%  fname_out = fs_label2label(fname,subj,'key1', value1,...); 
%
% Required Parameters:
%   fname:  full pathname of input label file
%
% Optional parameters:
%  'fname_out': output file name
%    If empty, will be derived from fname, with subj appended
%    {default = []}
%  'outdir': output directory
%    ignored if fname_out is full path
%    if empty, will use path of fname
%    {default = []}
%  'subj': FreeSurfer subject for whom output label will be made
%    {default = 'fsaverage'}
%  'source_subj': FreeSurfer subject from whom labels were made
%    {default = 'fsaverage'}
%  'subjdir': FreeSurfer subject root directory
%    {default = $FREESURFER_HOME/subjects}
%  'ico': icosahedron order number
%    used only if source_subj = 'fsaverage' or 'ico'
%    {default = 7}
%  'source_ico': icosahedron order number
%    {default = 7}
%  'verbose': [0|1] display mri_surf2surf output
%    {default = 1}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Created:  11/04/11 by Don Hagler
% Last Mod: 03/15/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin, 2), return; end;
parms = mmil_args2parms(varargin,{...
  'fname_in',fname,[],...
...
  'fname_out',[],[],...
  'outdir',[],[],...            
  'subj','fsaverage',[],...
  'source_subj','fsaverage',[],...
  'subjdir',getenv('SUBJECTS_DIR'),[],...
  'ico',7,[1:7],...
  'source_ico',7,[1:7],...
  'verbose',true,[false true],...
  'forceflag',false,[false true],...
...
  'hemi',[],{'lh','rh'},...
  'intype','label',[],...
  'outtype','label',[],...
});
fname_out = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check the input file name, extract the instem
[indir,inname,inext] = fileparts(parms.fname_in);
if ~strcmp(['.' parms.intype],inext)
  error('input fname %s missing correct file extension for intype %s',...
    parms.fname_in,parms.intype);
end;
pat = sprintf('(?<hemi>[lr]h).(?<instem>.+)');
n = regexp(inname,pat,'names');
if isempty(n)
  if isempty(parms.hemi)
    error('input file name %s missing hemi tag (lh or rh) at beginning of name and hemi parameter not specified',inname);
  else
    fprintf('%s: WARNING: input file name %s missing hemi tag (lh or rh) at beginning of name',mfilename,inname);
  end;
  instem = inname;
else
  instem = n.instem;
  parms.hemi = n.hemi;
end;
if isempty(parms.outdir)
  parms.outdir = indir;
end;

if isempty(parms.fname_out)
  % construct the output file name
  if (isempty(parms.fname_out))
    outstem = instem;
    if strcmp(indir,parms.outdir)
      if ismember(parms.subj,{'fsaverage','ico'})
        outstem = sprintf('%s-ico%d',outstem,parms.ico);
      else
        outstem = sprintf('%s-%s',outstem,parms.subj);
      end;
    end;
    parms.fname_out = sprintf('%s/%s.%s.%s',parms.outdir,parms.hemi,outstem,parms.outtype);
  end;
else
  if mmil_isrelative(parms.fname_out)
    parms.fname_out = sprintf('%s/%s',parms.outdir,parms.fname_out);
  else
    [parms.outdir,outname,outext] = fileparts(parms.fname_out);
  end;
end;
mmil_mkdir(parms.outdir);

if parms.forceflag || ~exist(parms.fname_out,'file')
    
  % construct the unix command
  cmd = sprintf('export SUBJECTS_DIR=''%s''\n',parms.subjdir); % cmd = sprintf('setenv SUBJECTS_DIR %s\n',parms.subjdir);
  cmd = sprintf('%smri_label2label',cmd);
  if ismember(parms.source_subj,{'fsaverage','ico'})
    cmd = sprintf('%s --srcsubject ico',cmd);
    cmd = sprintf('%s --srcicoorder %d',cmd,parms.source_ico);
  else
    cmd = sprintf('%s --srcsubject %s',cmd,parms.source_subj);
  end;
  if ismember(parms.subj,{'fsaverage','ico'})
    cmd = sprintf('%s --trgsubject ico',cmd);
    cmd = sprintf('%s --trgicoorder %d',cmd,parms.ico);
  else
    cmd = sprintf('%s --trgsubject %s',cmd,parms.subj);
  end;
  cmd = sprintf('%s --hemi %s',cmd,parms.hemi);
  cmd = sprintf('%s --srclabel %s',cmd,parms.fname_in);
  cmd = sprintf('%s --trglabel %s',cmd,parms.fname_out);
  cmd = sprintf('%s --regmethod surface',cmd);
  if parms.verbose
    fprintf('%s: cmd:\n%s\n',mfilename,cmd);
  end;
  % run the unix command
  [status,result]=unix(cmd);
  % check for success or failure
  if (status ||...
       ~isempty(findstr(result, 'could not open')) ||...
       ~exist(parms.fname_out,'file'))
    error('mri_label2label failed:\n%s\ncmd:\n%s',result,cmd);
  end;
  if parms.verbose, fprintf('%s',result); end;
end;

fname_out = parms.fname_out;

