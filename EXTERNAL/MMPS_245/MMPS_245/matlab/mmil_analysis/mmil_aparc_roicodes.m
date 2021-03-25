function roicodes = mmil_aparc_roicodes(varargin)
%function roicodes = mmil_aparc_roicodes([options])
%
% Usage:
%  errcode = mmil_aseg_roicodes('key1', value1,...);
%
% Optional Parameters:
%   'aparc_roilist': vector of aparc ROI codes to include
%     {default = [1001:1003,1005:1034,2001:2003,2005:2034]}
%   'exclude_roilist': vector of ROI codes to exclude
%     {default = []}
%
% Created:  08/29/17 by Don Hagler
% Last Mod: 08/29/17 by Don Hagler
%

%% todo: option for Desikan vs Destrieux parcellations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parse input arguments
parms = mmil_args2parms(varargin,{...
  'aparc_roilist',[1001:1003,1005:1034,2001:2003,2005:2034],[1,Inf],...
  'exclude_roilist',[],[1,Inf],...
});

% set ROI code list
roicodes = parms.aparc_roilist;
roicodes = setdiff(roicodes,parms.exclude_roilist);

