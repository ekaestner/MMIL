function roicodes = mmil_aseg_roicodes(varargin)
%function roicodes = mmil_aseg_roicodes([options])
%
% Usage:
%  errcode = mmil_aseg_roicodes('key1', value1,...);
%
% Optional Parameters:
%   'aseg_aparc_flag': [0|1] whether to use cortical parcellation ROIs
%     0: aseg only
%     1: aparc only   NOTE: for best results with aparc, use erode_flag=0
%     2: aparc+aseg
%     {default = 0}
%   'aseg_roilist': vector of aseg ROI codes to include
%     if aseg_aparc_flag = 0 or 2
%     {default = [1:28,40:60]}
%   'aparc_roilist': vector of aparc ROI codes to include
%     if aseg_aparc_flag = 1 or 2
%     {default = [1001:1003,1005:1034,2001:2003,2005:2034]}
%   'exclude_roilist': vector of ROI codes to exclude
%     {default = [1,3,6,9,19:23,25,27,40,42,45,48,55,56,57,59]}
%
% Created:  10/08/12 by Don Hagler
% Last Mod: 06/03/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parse input arguments
parms = mmil_args2parms(varargin,{...
  'aseg_aparc_flag',0,[0,1,2],...
  'aseg_roilist',[1:28,40:60],[1,Inf],...
  'aparc_roilist',[1001:1003,1005:1034,2001:2003,2005:2034],[1,Inf],...
  'exclude_roilist',[1,3,6,9,19:23,25,27,40,42,45,48,55,56,57,59],[1,Inf],...
});

% set ROI code list
switch parms.aseg_aparc_flag
  case 0
    roicodes = parms.aseg_roilist;
  case 1
    roicodes = parms.aparc_roilist;
  case 2
    roicodes = [parms.aseg_roilist,parms.aparc_roilist];
end;
roicodes = setdiff(roicodes,parms.exclude_roilist);

