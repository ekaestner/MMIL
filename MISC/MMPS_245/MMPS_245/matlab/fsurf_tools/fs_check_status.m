function [allexist,volexist,missing_files] = fs_check_status(subj,subjdir,version)
%function [allexist,volexist,missing_files] = fs_check_status(subj,[subjdir],[version])
%
% Required input:
%  subj is a string specifying the subject name
%
% Optional parameters:
%  subjdir - subjects directory (override SUBJECTS_DIR environment variable)
%    subjdir/subj should contain the freesurfer subject directory
%    {default = $SUBJECTS_DIR}
%  version - freesurfer version number (e.g. 305, 450, etc.)
%    {default = 450}
%
% Output:
%  allexist: returns 1 if all touchfiles exist, 0 otherwise
%  volexist: returns 1 if all volume touchfiles exist, 0 otherwise
%  missing_files: a cell array of the missing touch file names
%
% created:  12/22/06 by Don Hagler
% last mod: 11/02/11 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check parameters
if ~mmil_check_nargs(nargin,1), return; end;
allexist = 0;
volexist = 0;

if nargout==3
  missing_flag = 1;
  missing_files = [];
else
  missing_flag = 0;
end;

if ~exist('subjdir','var'), subjdir = []; end;
if isempty(subjdir)
  subjdir = getenv('SUBJECTS_DIR');
  if isempty(subjdir)
    error('SUBJECTS_DIR not defined as an environment variable');
  end;
end;

touchdir = sprintf('%s/%s/touch',subjdir,subj);
if ~exist(touchdir,'dir')
  fprintf('%s: WARNING: touch dir %s not found\n',mfilename,touchdir);
  return;
end

if ~exist('version','var') | isempty(version), version = 450; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if version >= 500 % e.g. version 5.1.0
  vol_steps = {...
    'motioncor' 'nuintensitycor' 'talairach' 'tal-check' ...
    'normalization' 'skullstrip' 'gcareg' 'canorm' 'careg' 'careginv' ...
    'rmneck' 'skull-lta' 'calabel' 'segstats'...
  };
  surf_steps = {...
    'normalization2' 'maskbfs' 'segmentation' 'fill' 'tessellate'...
    'smooth1' 'inflate1' 'qsphere' 'fix' 'white' 'smooth2' 'inflate2' ...
    'curvstats' 'sphere' 'surfreg' 'jacobian_white' ...
    'avgcurv' 'cortparc' 'pial' 'parcstats' 'cortparc2' 'parcstats2' ...
    'cortribbon' 'segstats' 'aparc2aseg' 'wmparc'...
    'label-exvivo-ec' 'ba-labels'...
  };
elseif version >= 400 % e.g. version 4.5.0
  vol_steps = {...
    'motioncor' 'nuintensitycor' 'talairach' 'tal-check' ...
    'normalization' 'skullstrip' 'gcareg' 'canorm' 'careg' 'careginv' ...
    'rmneck' 'skull-lta' 'calabel' 'segstats'...
  };
  surf_steps = {...
    'normalization2' 'maskbfs' 'segmentation' 'fill' 'tessellate'...
    'smooth1' 'inflate1' 'qsphere' 'fix' 'finalsurfs' 'smooth2' 'inflate2' ...
    'curvstats' 'sphere' 'surfreg' 'jacobian_white' ...
    'avgcurv' 'cortparc' 'parcstats' 'cortparc2' 'parcstats2' ...
  };
%    'cortribbon' 'aparc2aseg' 'wmparc'... % doesn't work on cluster
else % e.g. versions 3.0.2, 3.0.5
  vol_steps = {...
    'motioncor' 'nuintensitycor' 'talairach' ...
    'normalization' 'skullstrip' 'gcareg' 'canorm' 'careg' 'careginv' ...
    'rmneck' 'skull-lta' 'calabel' 'segstats'...
  };
  surf_steps = {...
    'normalization2' 'maskbfs' 'segmentation' 'fill' 'tessellate'...
    'smooth1' 'inflate1' 'qsphere' 'fix' 'finalsurfs' 'smooth2' 'inflate2' ...
    'sphere' 'surfreg' 'contrasurfreg' ...
    'avgcurv' 'cortparc' 'parcstats' 'cortparc2' 'parcstats2' ...
    'cortribbon' 'aparc2aseg' 'wmparc'...
  };
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[volexist,missing_files] = check_files(vol_steps,touchdir,missing_flag);
if ~volexist, return; end;
[allexist,missing_files] = check_files(surf_steps,touchdir,missing_flag);

return;

function [allexist,missing_files] = check_files(flist,touchdir,missing_flag)
  allexist = 1;
  missing_files = [];
  for f=1:length(flist)
    fname = [touchdir '/fs.finish.' flist{f} '.touch']; 
    if ~exist(fname,'file')
      allexist = 0;
      if missing_flag
        missing_files{end+1} = fname;
      else
        return;
      end;
    end;
  end;
return;
