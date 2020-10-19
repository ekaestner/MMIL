function [orient,oblique_flag] = fs_read_orient(fname,M)
%function [orient,oblique_flag] = fs_read_orient(fname,[M])
%
% Purpose: read header of mgh file and determine volume orientation
%
% Required Input:
%   fname: full path name of mgh/mgz file
%    If empty, expects M to be supplied
%
% Optional Input:
%   M: vox2ras matrix
%    {default = []}
%
% Output:
%   orient: three character string specifying volume orientation
%     e.g. 'RAS', 'LPI', 'PRI', etc.
%   oblique_flag: [0|1] whether slices are oblique
%
% Created:  02/14/10 by Don Hagler
% Last Mod: 05/03/10 by Don Hagler
%

if (~mmil_check_nargs(nargin,1)) return; end;
orient = 'XXX';
oblique_flag = 0;
smf = 10*eps;

if ~exist('M','var') || isempty(M)
  [vol,M] = fs_load_mgh(fname,[],[],1);
end;

M = M(1:3,1:3);
voxsize = sqrt(sum(abs(M).^2,1));

Rvec = [1 0 0]*M./voxsize;
Avec = [0 1 0]*M./voxsize;
Svec = [0 0 1]*M./voxsize;

[tmp,Rdim] = max(abs(Rvec));
[tmp,Adim] = max(abs(Avec));
[tmp,Sdim] = max(abs(Svec));

if Rvec(Rdim)>0
  orient(Rdim) = 'R';
else
  orient(Rdim) = 'L';
end;

if Avec(Adim)>0
  orient(Adim) = 'A';
else
  orient(Adim) = 'P';
end;

if Svec(Sdim)>0
  orient(Sdim) = 'S';
else
  orient(Sdim) = 'I';
end;

if nargout>1
  if length(find(abs(Rvec)>smf))>1 ||...
    length(find(abs(Avec)>smf))>1 ||...
    length(find(abs(Svec)>smf))>1
    oblique_flag = 1;
  end;
end;

