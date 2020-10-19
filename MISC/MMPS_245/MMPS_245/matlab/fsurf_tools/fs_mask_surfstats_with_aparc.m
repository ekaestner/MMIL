function masked_surfstats=fs_mask_surfstats_with_aparc(surfstats,subj,hemi,subjdir,mask_roinums)
%function masked_surfstats=fs_mask_surfstats_with_aparc(surfstats,subj,hemi,[subjdir],[mask_roinums])
% 
% Purpose: mask out portion of surface map based on aparc roi(s)
% 
% Required input:
%   surfstats: 2D matrix (nverts x frames) of surface stats
%   subj: freesurfer subject name
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
% Output:
%  masked_surfstats: input values masksed by aparc
%
% created:        02/12/07 Don Hagler
% last modified:  04/18/07 Don Hagler
%

masked_surfstats = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check parms

if nargin<3
  help(mfilename);
  return;
end;

if ~exist('subjdir','var'), subjdir = []; end;
if isempty(subjdir)
  subjdir = getenv('SUBJECTS_DIR');
  if isempty(subjdir)
    fprintf('%s: ERROR: SUBJECTS_DIR not defined as an environment variable\n',mfilename);
    return;
  end;
end;
if ~exist('mask_roinums','var') | isempty(mask_roinums), mask_roinums=[1,5]; end;

aparcname = sprintf('%s/%s/label/%s.aparc.annot',subjdir,subj,hemi);
if ~exist(aparcname,'file')
  fprintf('%s: ERROR: aparc file %s not found\n',...
    mfilename,aparcname);
  return;
end;      
[roinums,roilabels] = fs_read_annotation(aparcname);
if isempty(roinums)
  fprintf('%s: ERROR: unable to read annotation file %s\n',...
    mfilename,aparcname);
  return;
end;      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

roi = find(ismember(roinums,mask_roinums));
if size(surfstats,1) ~= length(roinums)
  fprintf('%s: ERROR: number of vertex data points (%d) does not match annotation file %s\n',...
    mfilename,size(surfstats,1),aparcname);
  return;
end;
masked_surfstats = surfstats;
masked_surfstats(roi,:) = 0;

