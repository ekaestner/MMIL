function RootDirs = MMIL_Set_Common_RootDirs(RootDir,fields)
%function RootDirs = MMIL_Set_Common_RootDirs(RootDir,[fields])
%
% Usage:
%  RootDirs = MMIL_Set_Common_RootDirs(RootDirs)
%
% Required Parameters:
%  RootDir: full path of root directory for all container types
%
% Optional Parameters:
%   fields: cell array of strings indicating which RootDirs fields
%    to set to RootDir
%    {default: {'orig','pc','qc','raw','proc','proc_dti','proc_bold',...
%               'fsurf','fsico','long','long_dti','orig_pet','raw_pet',...
%               'proc_pet','orig_meg','raw_meg','proc_meg','fiberdir'}}
%
% Created:  10/15/09 by Don Hagler
% Last Mod: 02/02/17 by Don Hagler
%

RootDirs = [];
if ~exist('fields','var') | isempty(fields)
  fields =  {'incoming','unpack','orig','pc','qc',...
             'raw','proc','proc_dti','proc_bold',...
             'fsurf','fsico','long','long_dti','orig_pet','raw_pet',...
             'proc_pet','orig_meg','raw_meg','proc_meg','fiberdir'};
end;

for i=1:length(fields)
  RootDirs = setfield(RootDirs,fields{i},RootDir);
end;

