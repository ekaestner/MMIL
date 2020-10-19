function [dimsmatch,orientmatch] = fs_vol_match(fname1,fname2)
%function [dimsmatch,orientmatch] = fs_vol_match(fname1,fname2)
%
% Purpose: read header of mgh files and check if volume dimensions match
%
% Required Input:
%   fname1: full path name of first mgh/mgz file
%   fname2: full path name of second mgh/mgz file
%
% Output:
%   dimsmatch: [0|1] whether first three dimensions match
%   orientmatch: [0|1] whether volume orientations match
%
% Created:  02/08/11 by Don Hagler
% Last Mod: 02/08/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,2), return; end;
dimsmatch = 0;
orientmatch = 0;

% check that input volumes share dimensions
[V,M1,mp,volsz1]=fs_load_mgh(fname1,[],[],1);
[V,M2,mp,volsz2]=fs_load_mgh(fname2,[],[],1);
diff_volsz = volsz1(1:3) - volsz2(1:3);
if all(diff_volsz==0), dimsmatch = 1; end;

if nargout>1
  orient1 = fs_read_orient([],M1);
  orient2 = fs_read_orient([],M2);
  if strcmp(orient1,orient2), orientmatch = 1; end;
end;

