function aparc_stats = fs_read_aparc_stats(subj,subjdir,aparc_infix);
%function aparc_stats = fs_read_aparc_stats(subj,[subjdir],[aparc_infix]);
%
% Required input:
%  subj: FreeSurferr recon subject name
%
% Optional input:
%  subjdir: subjects directory (override SUBJECTS_DIR environment variable)
%    subjdir/subj should contain the freesurfer subject directory
%    {default = $SUBJECTS_DIR}
%  aparc_infix: read {hemi}.{aparc_infix}.stats file from stats dir
%    {default = 'aparc'}
%
% Output:
%   aparc_stats: struct array containing:
%     roiname: name of ROI, with 'ctx-lh-' or 'ctx-rh-' preprended
%     roicode: numerical ROI code
%     hemi: hemisphere ('lh' or 'rh')
%     grayvol
%     surfarea
%     thickavg
%     thickstd
%     meancurv
%     gausscurv
%     foldind
%     curvind
%
% Created:  04/01/09 by Don Hagler
% Last Mod: 10/08/13 by Don Hagler
%

%% todo: roigroups for lobar analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aparc_stats = [];
if ~mmil_check_nargs(nargin,1), return; end;

hemilist = {'lh','rh'};

if ~exist('aparc_infix','var') | isempty(aparc_infix)
  aparc_infix = 'aparc';
end;

if ~exist('subjdir','var') | isempty(subjdir)
  subjdir = getenv('SUBJECTS_DIR');
  if isempty(subjdir)
    error('SUBJECTS_DIR not defined as environment variable');
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for h=1:length(hemilist)
  hemi = hemilist{h};
  tmp_stats = read_stats_file(subj,subjdir,aparc_infix,hemi);
  aparc_stats = [aparc_stats,tmp_stats];
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stats = read_stats_file(subj,subjdir,aparc_infix,hemi)
  stats = [];
  fname_aparc = sprintf('%s/%s/stats/%s.%s.stats',...
    subjdir,subj,hemi,aparc_infix);
  if ~exist(fname_aparc,'file')
    fprintf('%s: WARNING: aparc stats file %s not found\n',...
      mfilename,fname_aparc);
    return;
  end;
  % get aparc roi codes and names
  [all_roicodes,all_roinames] = fs_colorlut;
  switch hemi
    case 'lh'
      [tmp,ind]=intersect(all_roicodes,[1000:1300]);
    case 'rh'
      [tmp,ind]=intersect(all_roicodes,[2000:2300]);
  end;
  roicodes = all_roicodes(ind);
  roinames = all_roinames(ind);
  nrois = length(ind);

  % read aparc file
  if ~isempty(fname_aparc)
    try
      fid = fopen(fname_aparc);
      tmp_stats = textscan(fid,'%s %d %d %d %f %f %f %f %f %f\n',...
        'commentstyle','#');
      for i=1:length(tmp_stats{1})
        roiname = char(tmp_stats{1}{i});
        ind = find(~cellfun(@isempty,...
          regexp(roinames,['^ctx-' hemi '-' roiname '$'])));
        if isempty(ind)
          roicode = NaN;
        else
          roicode = roicodes(ind);
          roiname = roinames{ind};
        end;
        stats(i).roiname   = roiname;
        stats(i).roicode   = roicode;
        stats(i).hemi      = hemi;
        stats(i).surfarea  = double(tmp_stats{3}(i));
        stats(i).grayvol   = double(tmp_stats{4}(i));
        stats(i).thickavg  = double(tmp_stats{5}(i));
        stats(i).thickstd  = double(tmp_stats{6}(i));
        stats(i).meancurv  = double(tmp_stats{7}(i));
        stats(i).gausscurv = double(tmp_stats{8}(i));
        stats(i).foldind   = double(tmp_stats{9}(i));
        stats(i).curvind   = double(tmp_stats{10}(i));
      end;
    catch
      fprintf('%s: WARNING: failed to read %s aparc stats file\n',...
        mfilename,hemi);
    end;
  end;
return;

