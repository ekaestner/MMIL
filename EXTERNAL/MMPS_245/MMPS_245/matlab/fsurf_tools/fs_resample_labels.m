function label_fnames = fs_resample_labels(label_fnames,varargin)
%function label_fnames = fs_resample_labels(label_fnames,[options])
%
% Purpose: resample label files to another subject
%
% Required Input:
%   label_fnames: cell array of FreeSurfer label file names
%     standard naming convention is lh.<labelname>.label and rh.<labelname>.label
%
% Optional Input:
%   'labelnames': cell array of label names
%     If label_fnames is empty, will use files with labelnames in indir
%     {default = []}
%   'indir': input directory
%     If label_fnames is empty, will use label files found in indir
%     Also used if file names in label_fnames are not absolute path
%     {default = pwd}
%   'outdir': output directory
%     {default = pwd}
%   'source_subj': FreeSurfer subject from whom labels were made
%     {default = 'fsaverage'}
%   'subj': FreeSurfer subject for whom annotation will be made
%     {default = 'fsaverage'}
%   'source_ico': icosahedral order if 'source_subj' is 'ico' or 'fsaverage'
%     {default = 7}
%   'ico': icosahedral order if 'subj' is 'ico' or 'fsaverage'
%     {default = 7}
%   'subjdir': FreeSurfer subject root directory
%     {default = $FREESURFER_HOME/subjects}
%   'infix': extra string added to output file names after <labelname>
%     e.g. lh.<labelname>.<infix>.label
%     {default = []}
%   'hemilist': cell array of hemispheres for which to resample labels
%     {default = {'lh','rh'}}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Output:
%   label_fnames: cell array of label files found or created
%
% Created:  11/16/11 by Don Hagler
% Rcnt Mod: 11/16/11 by Don Hagler
% Last Mod: 09/15/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin,{...
  'label_fnames',label_fnames,[],...
  'labelnames',[],[],...
  'indir',pwd,[],...
  'outdir',pwd,[],...
  'source_subj','fsaverage',[],...
  'subj','fsaverage',[],...
  'source_ico',7,[1:7],...
  'ico',7,[1:7],...
  'subjdir',[],[],...
  'infix',[],[],...
  'hemilist',{'lh','rh'},{'lh','rh'},...
  'forceflag',false,[false true],...
...
  'verbose',true,[false true],...
});
label_fnames = [];
nhemi = length(parms.hemilist);

if ~exist(parms.indir,'dir')
  error('indir %s not found',parms.indir);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find label files
if ~isempty(parms.label_fnames) % check label_fnames
  if ~iscell(parms.label_fnames), parms.label_fnames = {parms.label_fnames}; end;
  for f=1:length(parms.label_fnames)
    if mmil_isrelative(parms.label_fnames{f})
      parms.label_fnames{f} = sprintf('%s/%s',parms.indir,parms.label_fnames{f});
    end;
    if ~exist(parms.label_fnames{f},'file')
      error('file %s not found',parms.label_fnames{f});
    end;
  end;
elseif ~isempty(parms.labelnames) % find label files in indir that match labelnames
  if ~iscell(parms.labelnames), parms.labelnames = {parms.labelnames}; end;
  f = 1;
  for i=1:length(parms.labelnames)
    for h=1:nhemi
      hemi = parms.hemilist{h};
      labelname = parms.labelnames{i};
      fname = sprintf('%s/%s.%s.label',parms.indir,hemi,labelname);
      if ~exist(fname,'file')
        fprintf('%s: WARNING: file %s not found',mfilename,fname);
      else
        parms.label_fnames{f} = fname;
        f = f + 1;
      end;
    end;  
  end;  
else % find all label files in indir
  f = 1;
  for h=1:nhemi
    hemi = parms.hemilist{h};
    flist = dir(sprintf('%s/%s.*.label',parms.indir,hemi));
    for k=1:length(flist)
      parms.label_fnames{f} = sprintf('%s/%s',parms.indir,flist(k).name);
      f = f + 1;
    end;
  end;
end;

% resample labels from source_subj to subj
if ~strcmp(parms.subj,parms.source_subj) ||...
   (ismember(parms.subj,{'fsaverage','ico'}) && parms.ico~=parms.source_ico)
  mmil_mkdir(parms.outdir); % will error if cannot create
  resample_flag = 1;
else
  resample_flag = 0;
end;

for f=1:length(parms.label_fnames)
  fname_in = parms.label_fnames{f};
  n = regexp(fname_in,'(?<hemi>[lr]h)\.(?<name>.+)\.label$','names');
  if isempty(n)
    fprintf('%s: WARNING: unexpected naming convention for %s\n',...
      mfilename,fname_in);
    [tpath,tstem,text] = fileparts(fname_in);
    fname_out = [parms.outdir '/' tstem text];
    if length(parms.hemilist)==1
      hemi = parms.hemilist{1};
    else
      hemi = [];
    end;
  else
    hemi = n.hemi;
    labelname = n.name;
    if ~isempty(parms.infix)
      fname_out = sprintf('%s/%s.%s.%s.label',parms.outdir,hemi,labelname,parms.infix);
    else
      fname_out = sprintf('%s/%s.%s.label',parms.outdir,hemi,labelname);
    end;
  end;
  if isempty(hemi) || ~ismember(hemi,parms.hemilist)
    fprintf('%s: WARNING: skipping file %s...\n',mfilename,fname_in);
    continue;
  end;
  if resample_flag
    if strcmp(fname_in,fname_out)
      fprintf('%s: WARNING: input and output file names are identical for %s\n',...
        mfilename,fname_in);
      continue;
    else
      % resample label from source_subj to subj
      tags = {'subj','source_subj','subjdir','ico','source_ico','verbose','forceflag'};
      args = mmil_parms2args(parms,tags);
      if ~exist(fname_out,'file') || parms.forceflag
        fname_out = fs_label2label(fname_in,'fname_out',fname_out,'hemi',hemi,args{:});
      end;
    end;
    label_fnames{end+1} = fname_out;
  else
    label_fnames{end+1} = fname_in;
  end;
end;

