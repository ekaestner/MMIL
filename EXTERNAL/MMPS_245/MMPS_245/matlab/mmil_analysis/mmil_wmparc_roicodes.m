function roicodes = mmil_wmparc_roicodes(varargin)
%function roicodes = mmil_wmparc_roicodes([options])
%
% Usage:
%  errcode = mmil_wmparc_roicodes('key1', value1,...);
%
% Optional Parameters:
%   'wmparc_aparc_flag': [0|1] whether to use cortical parcellation ROIs
%     0: wmparc only
%     1: wmparc and aparc
%     {default = 0}
%   'wmparc_roilist': vector of aseg ROI codes to include
%     {default = [1:28,40:60]}
%   'aparc_roilist': vector of aparc ROI codes to include
%     ignored if wmparc_aparc_flag = 0
%     {default = [1001:1003,1005:1034,2001:2003,2005:2034]}
%   'exclude_roilist': vector of ROI codes to exclude
%     {default = [1,3,6,9,19:23,25,27,40,42,45,48,55,56,57,59]}
%
% Created:  06/03/12 by Don Hagler
% Last Mod: 06/03/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parse input arguments
parms = mmil_args2parms(varargin,{...
  'wmparc_aparc_flag',false,[false,true],...
  'wmparc_roilist',[3001:3003,3005:3034,4001:4003,4005:4034],[1,Inf],...
  'aparc_roilist',[1001:1003,1005:1034,2001:2003,2005:2034],[1,Inf],...
  'exclude_roilist',[],[1,Inf],...
});

% set ROI code list
if ~parms.wmparc_aparc_flag
  roicodes = parms.wmparc_roilist;
else
  roicodes = [parms.wmparc_roilist,parms.aparc_roilist];
end;
roicodes = setdiff(roicodes,parms.exclude_roilist);

