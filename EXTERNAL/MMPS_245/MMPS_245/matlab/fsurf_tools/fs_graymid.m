function fs_graymid(subj,subjdir,forceflag)
%function fs_graymid(subj,[subjdir],[forceflag])
%
% Purpose: generate a surface in the middle of gray-matter
%
% Required Input:
%  subj: subject name
%
% Optional Input:
%  subjdir: subjects directory
%     {default: SUBJECTS_DIR}
%  forceflag: [0|1] whether to force overwrite of existing graymid files
%     {default: 0}
%
% Created:   10/06/09 by Don Hagler
% Last Mod:  06/18/14 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~mmil_check_nargs(nargin,1)) return; end;
if ~exist('forceflag','var') | isempty(forceflag), forceflag = 0; end;
if ~exist('subjdir','var') | isempty(subjdir)
  subjdir = getenv('SUBJECTS_DIR');
  if isempty(subjdir), error('subdir not set and missing SUBJECTS_DIR'); end;
end;
hemilist = {'lh','rh'};
surfdir = [subjdir '/' subj '/surf'];

for h=1:length(hemilist)
  hemi = hemilist{h};
  fname_out = [surfdir '/' hemi '.graymid'];
  if ~exist(fname_out,'file') | forceflag
    cmd = ['cd ' surfdir '; mris_expand -thickness ' hemi '.white 0.5 ' hemi '.graymid'];
    fprintf('%s: cmd = %s\n',mfilename,cmd);
    [s,r] = unix(cmd);
    if s
      error('cmd %s failed:\n%s',cmd,r);
    end;
    fprintf('%s\n',r);
  end;
end;

