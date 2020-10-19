function [vol,qmat,bvals,ind_censored_frames] =...
    dti_censor_frames(fname_censor,vol,qmat,bvals)
%function [vol,qmat,bvals,ind_censored_frames] =...
%    dti_censor_frames(fname_censor,[vol],[qmat],[bvals])
%
% Purpose: remove frames from diffusion data volume, qmat, and bvals
%   according to censor.txt file
%
% Required Parameters:
%   fname_censor: full path name of text file containing frame numbers (1-based)
%     should be arranged in single column
%
% Optional Parameters:
%   vol: diffusion MRI volume
%   qmat: matrix of diffusion direction vectors (ndirs x 3)
%   bvals: vector of b-values (can be single value if same for all)
%
% Note: vol or qmat can be empty, but not both
%       bvals can be empty
%
% Created:  11/23/10 by Don Hagler
% Last Mod: 10/27/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin, 1), return; end;

ind_censored_frames = [];
if ~exist('vol','var'), vol = []; end;
if ~exist('qmat','var'), qmat = []; end;
if ~exist('bvals','var'), bvals = []; end;
if ~isempty(fname_censor) & exist(fname_censor,'file')
  if isempty(vol) && isempty(qmat)
    error('vol and qmat cannot both be empty');
  elseif isempty(vol)
    nf = size(qmat,1);
  else
    nf = size(vol,4);
  end;
  ind_censored_frames = textread(fname_censor);
  ind_censored_frames = intersect([1:nf],ind_censored_frames);
  fprintf('%s: %d of %d censored frames from %s:\n',...
    mfilename,length(ind_censored_frames),nf,fname_censor);
  disp(ind_censored_frames);
  ind_keep = setdiff([1:nf],ind_censored_frames);
  if ~isempty(vol)
    vol = vol(:,:,:,ind_keep);
  end;
  if ~isempty(qmat)
    qmat = qmat(ind_keep,:);
  end;
  if ~isempty(bvals) && length(bvals)==nf
    bvals = bvals(ind_keep);
  end;
end;


