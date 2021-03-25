function ABCD_Combine_Import(ProjID,varargin)
%function ABCD_Combine_Import(ProjID,[options])
%
% Combine import_info from mmilrec14 and mmilrec18
%
% Required Input:
%   ProjID: project ID string
%
%
% Created:  03/27/17 by Feng Xue
% Last Mod: 03/27/17 by Feng Xue
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
abcd_combine_import;
