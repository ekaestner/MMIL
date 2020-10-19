function fs_mask_surfmgh_with_aparc(subj,hemi,infile,outfile,subjdir,mask_roinums)
%function fs_mask_surfmgh_with_aparc(subj,hemi,infile,outfile,[subjdir],[mask_roinums])
% 
% Purpose: mask out portion of surface map based on aparc roi(s)
% 
% Required input:
%   subj: freesurfer subject name
%   infile: full/relative path of input surface map (must be mgh format)
%   outfile: full/relative path of output surface map (must be mgh format)
%   hemi: cortical hemisphere ('lh' or 'rh')
%
% Optional parameters:
%   subjdir: subjects directory (override SUBJECTS_DIR environment variable)
%    subjdir/subj should contain the freesurfer subject directory
%    {default = $SUBJECTS_DIR}
%   mask_roinums: vector of numbers between 1 and 35 corresponding to aparc ROIs
%    e.g. 1=unknown, 5=corpuscallosum (see subjdir/subj/stats/lh.aparc.stats)
%    {default = [1,5]}
%
% created:        02/12/07 Don Hagler
% last modified:  09/18/07 Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check parms

if nargin<4
  help(mfilename);
  return;
end;

if ~exist('subjdir','var'), subjdir = []; end;
if isempty(subjdir)
  subjdir = getenv('SUBJECTS_DIR');
  if isempty(subjdir)
    error('SUBJECTS_DIR not defined as an environment variable');
  end;
end;
if ~exist('mask_roinums','var') | isempty(mask_roinums), mask_roinums=[1,5]; end;

if ~exist(infile,'file')
  error('input file %s not found',infile);
end;      

aparcname = sprintf('%s/%s/label/%s.aparc.annot',subjdir,subj,hemi);
if ~exist(aparcname,'file')
  error('aparc file %s not found',aparcname);
end;      
[roinums,roilabels] = fs_read_annotation(aparcname);
if isempty(roinums)
  error('unable to read annotation file %s',aparcname);
  return;
end;      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

roi = find(ismember(roinums,mask_roinums));
[vol,M] = fs_load_mgh(infile);
if isempty(vol)
  error('input volume is empty');
end;

if length(size(vol))==2
  nverts = prod(size(vol));
  vol = reshape(vol,[nverts,1,1,1]);
else
  nverts = size(vol,1);
end;
if nverts ~= length(roinums)
  error('number of vertex data points in %s (%d) does not match annotation file %s (%d)',...
    infile,nverts,aparcname,length(roinums));
end;
vol(roi,:,:,:) = 0;
fs_save_mgh(vol,outfile,M);

