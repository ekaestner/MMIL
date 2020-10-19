function [label_fnames,ctab] = fs_annot2labels(fname,varargin)
%function [label_fnames,ctab] = fs_annot2labels(fname,[options])
%
% Purpose: convert an annotation file into a collection of label files
%
% Required Input:
%   fname: full or relative path of annotation file
%
% Optional Input:
%   'outdir': output directory
%     {default = pwd}
%   'subj': FreeSurfer subject name included in label file
%     {default = 'fsaverage'}
%   'hemi': hemisphere string ('lh' or 'rh')
%     not required if fname has form {hemi}.{name}.annot
%     {default = []}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Output:
%   label_fnames: cell array of filenames for left and right hemispheres
%   ctab: color table struct
%     ctab.numEntries = number of Entries
%     ctab.orig_tab = name of original ct
%     ctab.struct_names = list of structure names (e.g. central sulcus and so on)
%     ctab.table = n x 5 matrix. 1st column is r, 2nd column is g, 3rd column
%       is b, 4th column is flag, 5th column is resultant integer values
%     calculated from r + g*2^8 + b*2^16 + flag*2^24. flag expected to be all 0
%
% Created:  11/29/12 by Don Hagler
% Last Mod: 09/11/15 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
label_fnames = []; ctab = [];

parms = mmil_args2parms(varargin,{...
  'outdir',pwd,[],...
  'subj','fsaverage',[],...
  'hemi',[],[],...
  'forceflag',false,[false true],...
...
  'verbose',true,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(fname,'file'), error('files %s not found',fname); end;

% determine hemi using regexp
[tpath,tstem,text] = fileparts(fname);
n = regexp(tstem,'(?<hemi>[lr]h)\.\w','names');
if isempty(n)
  fprintf('%s: WARNING: fname %s has unexpected pattern...\n',...
    mfilename,fname);
else
  parms.hemi = n.hemi;
end;
if isempty(parms.hemi)
  error('hemi not specified');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read annot file
[roinums,roilabels,ctab] = fs_read_annotation(fname);

% save each label file
nlabels = length(roilabels);
j = 1;
for i=1:nlabels
  roiname = roilabels{i};
  if strcmp(roiname,'unknown'), continue; end;
  fname_out = sprintf('%s/%s.%s.label',parms.outdir,parms.hemi,roiname);
  if ~exist(fname_out,'file') || parms.forceflag
    mmil_mkdir(parms.outdir);
    vertices = find(roinums == i);
    if ~isempty(vertices)
      fs_write_label(vertices,fname_out,parms.subj);
    end;
  end;
  if exist(fname_out,'file')
    label_fnames{j} = fname_out;
    j = j + 1;
  end;
end;

