function [M,volsz,orient,mr_parms] = fs_read_header(fname)
%function [M,volsz,orient,mr_parms] = fs_read_header(fname)
%
% Purpose: read header of mgh file and determine volume orientation
%
% Required Input:
%   fname: full path name of mgh/mgz file
%
% Output:
%   M: vox2ras matrix
%   volsz: volume dimensions [rows,columns,slices,frames]
%   orient: three character string specifying volume orientation
%     e.g. 'RAS', 'LPI', 'PRI', etc.
%   mr_parms: [tr flipangle te ti fov]
%
% Created:  04/23/10 by Don Hagler
% Last Mod: 04/23/10 by Don Hagler
%

if (~mmil_check_nargs(nargin,1)) return; end;

[vol,M,mr_parms,volsz] = fs_load_mgh(fname,[],[],1);
orient = fs_read_orient([],M);

