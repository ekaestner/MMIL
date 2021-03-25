function [annot_fnames,label_fnames] = fs_labels2annot(labeldir,varargin)
%function annot_fnames = fs_labels2annot(labeldir,[options])
%
% Purpose: resample a directory of labels into an annotation file
%
% Required Input:
%   labeldir: full path of directory containing FreeSurfer label files
%     with naming convention being lh.<labelname>.label and rh.<labelname>.label
%
% Optional Input:
%   'fname_ctab': name of color table file with names matching <labelname>
%      if empty, will use $FREESURFER_HOME/FreeSurferColorLUT.txt
%     {default = []}
%   'outdir': output directory
%     {default = pwd}
%   'source_subj': FreeSurfer subject from whom labels were made
%     {default = 'fsaverage'}
%   'subj': FreeSurfer subject for whom annotation will be made
%     {default = 'fsaverage'}
%   'subjdir': FreeSurfer subject root directory
%     {default = $FREESURFER_HOME/subjects}
%   'annotname': output file stem
%     {default = 'fparc'}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Output:
%   annot_fnames: cell array of filenames for left and right hemispheres
%   label_fnames: cell array of label files found or created
%
% Created:  11/03/11 by Don Hagler
% Last Mod: 11/29/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
annot_fnames = {};
label_fnames = {};
parms = mmil_args2parms(varargin,{...
  'indir',labeldir,[],...
  'outdir',pwd,[],...
  'source_subj','fsaverage',[],...
  'fname_ctab',[],[],...
  'subj','fsaverage',[],...
  'subjdir',[],[],...
  'annotname','fparc',[],...
  'forceflag',false,[false true],...
...
  'hemilist',{'lh','rh'},{'lh','rh'},...
  'ico',7,[1:7],...
  'source_ico',7,[1:7],...
  'verbose',true,[false true],...
  'cleanup_flag',true,[false true],...
});

nhemi = length(parms.hemilist);
Ns = [42 162 642 2562 10242 40962 163842];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if mmil_isrelative(parms.outdir)
  parms.outdir = [pwd '/' parms.outdir];
end;

if isempty(parms.subjdir)
  parms.subjdir = sprintf('%s/subjects',getenv('FREESURFER_HOME'));
end;

if isempty(parms.fname_ctab)
  parms.fname_ctab = sprintf('%s/FreeSurferColorLUT.txt',...
                      getenv('FREESURFER_HOME'));
end;

% read ctab file
[roicodes,labelnames,rgbv] = fs_colorlut(parms.fname_ctab);
nroi = length(labelnames);

% construct ctab struct
ctab = struct();
ctab.numEntries = nroi;
ctab.orig_tab = parms.fname_ctab;
ctab.struct_names = labelnames;
ctab.table = zeros(ctab.numEntries,5);
ctab.table(:,1:4) = rgbv;
for i=1:ctab.numEntries
  ctab.table(i,5) =  ctab.table(i,1) + ...
                     ctab.table(i,2)*2^8 + ...
                     ctab.table(i,3)*2^16 + ...
                     ctab.table(i,4)*2^24;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check that label files exist for each labelname

label_fnames = {};
for h=1:nhemi
  hemi = parms.hemilist{h};
  for r=1:nroi
    labelname = labelnames{r};
    if strcmp(labelname,'unknown'), continue; end;
    fname_in = sprintf('%s/%s.%s.label',...
      parms.indir,hemi,labelname);
    if ~exist(fname_in,'file') && parms.verbose
      fprintf('%s: WARNING: file %s not found\n',mfilename,fname_in);
      continue;
    end;
    label_fnames{end+1} = fname_in;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% resample labels from source_subj to subj

tmpdir = sprintf('%s/tmp_labels_subj',parms.outdir);
tmp_parms = parms;
tmp_parms.outdir = tmpdir;
tags = {'subj','source_subj','subjdir','outdir',...
  'ico','source_ico','verbose','forceflag'};
args = mmil_parms2args(tmp_parms,tags);
label_fnames = fs_resample_labels(label_fnames,args{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create annot file from set of labels

mmil_mkdir(parms.outdir); % will error if cannot create
for h=1:nhemi
  hemi = parms.hemilist{h};

  % find number of vertices for entire hemisphere
  if ismember(parms.subj,{'fsaverage','ico'})
    nverts = Ns(parms.ico);
  else
    fname_surf = sprintf('%s/%s/surf/%s.white',parms.subjdir,parms.subj,hemi);
    nverts = fs_read_surf_nverts(fname_surf);
  end;

  fname_out = sprintf('%s/%s.%s.annot',parms.outdir,hemi,parms.annotname);
  if parms.forceflag || ~exist(fname_out,'file')
    colorcodes = zeros(nverts,1);
    % load labels
    for f=1:length(label_fnames)
      fname_in = label_fnames{f};
      n = regexp(fname_in,'(?<hemi>[lr]h)\.(?<name>.+)\.label$','names');
      if isempty(n)
        fprintf('%s: WARNING: skipping label %s with nonstandard naming convention\n',...
          mfilename,fname_in);
        continue;
      end;
      if ~strcmp(hemi,n.hemi), continue; end;
      labelname = n.name;
      ind_label = find(strcmp(labelnames,labelname));
      v = fs_read_label(fname_in);
      colorcodes(v) = ctab.table(ind_label,5);
    end;
    ind_unknown = find(strcmp(labelnames,'unknown'));
    if ~isempty(ind_unknown)
      v_unknown = setdiff(1:nverts,find(colorcodes));
      colorcodes(v_unknown) = ctab.table(ind_unknown,5);
    end;
    vertices = [0:nverts-1];
    fs_write_annotation(fname_out,vertices,colorcodes,ctab);
  end;
  annot_fnames{h} = fname_out;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% remove temporary files
if parms.cleanup_flag && exist(tmpdir,'dir')
  cmd = sprintf('rm -r %s',tmpdir);
  [s,r] = unix(cmd);
  if s, error('cmd %s failed:\n%s',cmd,r); end;
end;
