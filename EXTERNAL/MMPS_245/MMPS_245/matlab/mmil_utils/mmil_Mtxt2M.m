function M=mmil_Mtxt2M(fname_reg,varargin)
%function M=mmil_Mtxt2M(fname)
%
% Purpose: Converts 4*4 transformation matrix from reg txt
%   to M (similar to rbreg)
%  first and second rows swapped and first and second columns swapped
%  compared to output of rbreg_vol2vol_icorr_mask
%
% Required Input:
%   fname_reg: registration matrix text file from mmil_reg output
%
% Optional Input:
%   Mreg: Provide Mreg matrix to convert instead of fname (fname should be [])
%
% Output:
%   M: Transformation matrix in format similar to rbreg  
%   
%
% Created:  02/02/12 by Vijay Venkatraman
% Last Mod: 02/02/12 by Vijay Venkatraman
%
if (~mmil_check_nargs(nargin,1)) return; end;
parms = mmil_args2parms(varargin, {...
  'Mreg',[],[],...
});

if exist(fname_reg,'file')
  M_tmp= textread(fname_reg);
elseif ~isempty(parms.Mreg)
  M_tmp= parms.Mreg;
else
  error('%s:missing Mreg or fname_reg',mfilename);
end;

M = M_tmp;
M(1,:) = M_tmp(2,:);
M(2,:) = M_tmp(1,:);
M_tmp = M ;
M(:,1) = M_tmp(:,2);
M(:,2) = M_tmp(:,1);

return;

