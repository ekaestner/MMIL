function fs_calc_area_ratio(subj,varargin);
%function fs_calc_area_ratio(subj,[options]);
%
% Usage:
%  fs_calc_area_ratio(subj,'key1', value1,...);
%
% Required input:
%  subj is a string specifying the subject name
%
% Optional parameters:
%  'surfname1' - surface name 1
%    {default: sphere}
%  'surfname2' - surface name 2
%    {default: sphere.reg}
%  'hemi' - should be either 'lh' for left hemisphere or 'rh' for right hemi
%    {default = both}
%  'outstem'  - output file stem (omit extension, hemi)
%    {default = 'surfname1_surfname2_arearatio'}
%  'outtype' - output file type ('mgh' or 'w')
%    {default: 'mgh'}
%  subjdir - subjects directory (override SUBJECTS_DIR environment variable)
%    subjdir/subj should contain the freesurfer subject directory
%    {default = $SUBJECTS_DIR}
%  'smooth' - number of smoothing steps
%    {default: 100}
%  'logflag' - [0|1] toggle output log of ratio
%    {default: 1}
%  'overwrite_flag' - [0|1] toggle overwrite existing output files
%    {default: 0}
%
% created:        11/17/06 Don Hagler
% last modified:  08/20/07 Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set defaults
DEFAULT_SURFNAME1 = 'sphere';
DEFAULT_SURFNAME2 = 'sphere.reg';
DEFAULT_OUTTYPE = 'mgh';
DEFAULT_SMOOTH = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse options

try
  options = varargin;
  for index = 1:length(options)
      if iscell(options{index}) & ~iscell(options{index}{1}), options{index} = { options{index} }; end;
  end;
  if ~isempty( varargin ), g=struct(options{:}); 
  else g = []; end;
catch
  error('calling convention {''key'', value, ... } error');
end;

if ~isfield(g,'hemi'), g.hemi = []; end;
if ~isfield(g,'outstem'), g.outstem = []; end;
if ~isfield(g,'outtype'), g.outtype = DEFAULT_OUTTYPE; end;
if ~isfield(g,'subjdir'), g.subjdir = []; end;
if ~isfield(g,'surfname1'), g.surfname1 = DEFAULT_SURFNAME1; end;
if ~isfield(g,'surfname2'), g.surfname2 = DEFAULT_SURFNAME2; end;
if ~isfield(g,'smooth'), g.smooth = DEFAULT_SMOOTH; end;
if ~isfield(g,'logflag'), g.logflag = 1; end;
if ~isfield(g,'overwrite_flag'), g.overwrite_flag = 0; end;

gfields = fieldnames(g);
for index=1:length(gfields)
   switch gfields{index}
   case {'hemi' 'outstem' 'outtype' ...
         'subjdir' 'surfname1' 'surfname2'...
         'smooth' 'logflag' 'overwrite_flag'},;
   otherwise, error([mfilename ': unrecognized option: ''' gfields{index} '''' ]);
   end;
end;

% get rid of options struct
hemi = g.hemi;
outstem = g.outstem;
outtype = g.outtype;
subjdir = g.subjdir;
surfname1 = g.surfname1;
surfname2 = g.surfname2;
smooth = g.smooth;
logflag = g.logflag;
overwrite_flag = g.overwrite_flag;
clear g;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check parameters
if nargin<1, help(mfilename); return; end;  

if ~exist('subjdir','var'), subjdir = []; end;
if isempty(subjdir)
  subjdir = getenv('SUBJECTS_DIR');
  if isempty(subjdir)
    fprintf('%s: SUBJECTS_DIR not defined as an environment variable... quitting\n',mfilename);
    return;
  end;
end;

if isempty(hemi)
  hemilist = {'lh','rh'};
else
  if ~ismember(hemi,{'lh','rh'})
    fprintf('%s: hemi must be lh or rh (is %s)\n',mfilename,hemi);
    return;
  end;
  hemilist = {hemi};
end;

for h=1:length(hemilist)
  hemi = hemilist{h};
  surffile = sprintf('%s/%s/surf/%s.%s',subjdir,subj,hemi,surfname1);
  if ~exist(surffile,'file')
    fprintf('%s: surface file %s not found... quitting\n',mfilename,surffile);
    return;
  end
  surffile = sprintf('%s/%s/surf/%s.%s',subjdir,subj,hemi,surfname2);
  if ~exist(surffile,'file')
    fprintf('%s: surface file %s not found... quitting\n',mfilename,surffile);
    return;
  end
end;
setenv('SUBJECTS_DIR',subjdir);

if isempty(outstem)
  outstem = sprintf('%s_%s_arearatio',surfname1,surfname2);
end;

if ~ismember(outtype,{'mgh','w'})
  fprintf('%s: outtype must be mgh or w (is %s)\n',mfilename,outtype);
  return;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for h=1:length(hemilist)
  outfile = sprintf('%s/%s/surf/%s-%s.%s',subjdir,subj,outstem,hemi,outtype);
  if ~exist(outfile,'file') | overwrite_flag
    fprintf('%s: generating %s...\n',mfilename,outfile);
    hemi = hemilist{h};
    surffile = sprintf('%s/%s/surf/%s.%s',subjdir,subj,hemi,surfname1);
    surf1 = fs_read_surf(surffile);
    surf1 = fs_find_neighbors(surf1);
    surf1 = fs_calc_triarea(surf1);

    surffile = sprintf('%s/%s/surf/%s.%s',subjdir,subj,hemi,surfname2);
    surf2 = fs_read_surf(surffile);
    surf2 = fs_find_neighbors(surf2);
    surf2 = fs_calc_triarea(surf2);

    area1 = surf1.vertex_area;
    area2 = surf2.vertex_area;

    area1_sm = fs_smooth(surf1,area1,smooth);
    area2_sm = fs_smooth(surf2,area2,smooth);

    tmp_ind = find(area1_sm==0);
    area1_sm(tmp_ind)=1;
    area_ratio = area2_sm./area1_sm;
    area_ratio(tmp_ind)=0;
    if logflag
      area_ratio(find(area_ratio==0))=1;
      area_ratio = log(area_ratio);
    end;

    if strcmp(outtype,'mgh')
      fs_save_mgh(area_ratio,outfile,eye(4));
    else
      v = find(area_ratio);
      w = area_ratio(v);
      fs_write_wfile(outfile,w,v);
    end;
  end;
end;

return;

