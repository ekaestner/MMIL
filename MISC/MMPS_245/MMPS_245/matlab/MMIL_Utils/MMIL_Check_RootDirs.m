function RootDirs = MMIL_Check_RootDirs(RootDirs,required_rootdirs,...
                                                 ntry,pdur)
%function RootDirs = MMIL_Check_RootDirs(RootDirs,[required_rootdirs],...
%                                                 [ntry],[pdur])
%
% Usage:
%  RootDirs = MMIL_Check_RootDirs(RootDirs)
%
% Required Parameters:
%  RootDirs: struct that may contain the following fields:
%     batch, incoming, unpack, orig,
%     raw, proc, proc_dti, proc_bold,
%     fsurf, fsico, long,
%     orig_meg, raw_meg, proc_meg
%     orig_pet, raw_pet,  proc_pet
%   these specify the full paths of root directories containing data containers
%     if not found, they will be given empty values
%   batch specifies where batchdirs will go
%     if unspecified, will go in /home/<user>/batchdirs
%
% Optional Parameters:
%   required_rootdirs: cell array of strings indicating which 
%     RootDirs fields must not be empty
%    {default = []}
%   ntry: number of attempts to check for existence of each RootDir
%    {defualt = 10}
%   pdur: duration (sec) of pause between attempts to check for existence
%    {default = 0.5}
%
% Created:  03/31/09 by Don Hagler
% Last Mod: 02/02/17 by Don Hagler
%

if ~exist('required_rootdirs','var'), required_rootdirs = []; end;
if ~exist('ntry','var') || isempty(ntry), ntry = 10; end;
if ~exist('pdur','var') || isempty(pdur), pdur = 0.5; end;

ContainerTypes = {'home','batch','incoming','unpack','orig','pc','qc',...
                  'raw','proc','proc_dti','proc_bold',...
                  'raw_asl','proc_asl',...
                  'fsurf','fsico','long','long_dti',...
                  'orig_meg','raw_meg','proc_meg',...
                  'orig_pet','orig_pet_pib','raw_pet','proc_pet',...
                  'fiberdir'};
                  
for i=1:length(ContainerTypes)
  if ~isfield(RootDirs,ContainerTypes{i})
    RootDirs = setfield(RootDirs,ContainerTypes{i},[]);
  end;
  RootDir = getfield(RootDirs,ContainerTypes{i});
  if ~isempty(RootDir)
    exist_flag = 0;
    for n=1:ntry
      if exist(RootDir,'dir')
        exist_flag = 1;
        break;
      else
        fprintf('%s: WARNING: failed to find rootdir %s, waiting %0.1f sec to check again...\n',...
          mfilename,RootDir,pdur);
        pause(pdur);
      end;
    end;
    if ~exist_flag
      error('RootDir %s not found',RootDir);
    end;
  end;
end;

if isempty(RootDirs.home)
  RootDirs.home = getenv('HOME');
end;
if isempty(RootDirs.batch)
  RootDirs.batch = [RootDirs.home '/batchdirs'];
end;  

if ~isempty(required_rootdirs)
  if ~iscell(required_rootdirs), required_rootdirs = {required_rootdirs}; end;
  for f=1:length(required_rootdirs)
    fieldname = required_rootdirs{f};
    if ~isfield(RootDirs,fieldname) || isempty(RootDirs.(fieldname))
      % if proc_dti or proc_bold are missing, set to proc
      if ismember(fieldname,{'proc_dti','proc_bold'})
        RootDirs.(fieldname) = getfield(RootDirs,'proc');
        if ~isempty(RootDirs.(fieldname)), continue; end;
      end;
      error('RootDirs is missing field %s',fieldname);
    end;
  end;  
end;

