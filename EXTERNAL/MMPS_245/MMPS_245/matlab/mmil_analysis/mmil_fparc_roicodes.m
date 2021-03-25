function roicodes = mmil_fparc_roicodes(varargin)
%function roicodes = mmil_fparc_roicodes([options])
%
% Usage:
%  errcode = mmil_fparc_roicodes('key1', value1,...);
%
% Optional Parameters:
%   'fparc_roilist': vector of fparc ROI codes to include
%     {default = [21101:21300,22101:22300]}
%   'exclude_roilist': vector of ROI codes to exclude
%     {default = []}
%
% Created:  08/29/17 by Don Hagler
% Last Mod: 08/29/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parse input arguments
parms = mmil_args2parms(varargin,{...
  'fparc_roilist',[21101:21300,22101:22300],[1,Inf],...
  'exclude_roilist',[],[1,Inf],...
});

% set ROI code list
roicodes = parms.fparc_roilist;
roicodes = setdiff(roicodes,parms.exclude_roilist);

