function MMIL_Check_Usage(ContainerInDirs)
%function MMIL_Check_Usage(ContainerInDirs)
%
% Inputs:
% ContainerInDirs: string or cell array of container input directories to check hard drive space on
%
%
% Created 07/09/2009 by Alain Koyama

maxusage = 95; % percent usage at which a warning message will be outputted

if ~exist('ContainerInDirs','var') | isempty(ContainerInDirs)
  fprintf('%s: Warning - No input container(s) specified\n',mfilename);
  return; 
end;

if ischar(ContainerInDirs) % if there is a single dir inputted as a string, convert to cell
  ContainerInDirs = {ContainerInDirs};
end

for i=1:length(ContainerInDirs) % loop through all output dirs to be checked
  ContainerInDir = ContainerInDirs{i};
  hd_space = char(regexp(ContainerInDir,'^/space/.+/\d{1,2}/','match')); 
  hd_home = char(regexp(ContainerInDir,'^/home/.+','match'));
  if ~isempty(hd_space) % if it's a /space partition
    partition = hd_space;
  elseif ~isempty(hd_home) % if it's the /home partition
    partition = hd_home;
  else
    fprintf('%s: Warning - unable to parse the path %s for a valid partition\n',...
      mfilename,ContainerInDir);
    continue;
  end    
  [status,result] = unix(sprintf('df -h %s',partition));
  if status
    fprintf('%s: Warning - unable to check disk space usage on %s: %s\n',...
      mfilename,partition,result);
    continue;
  end
  usage = str2num(char(regexp(result,'\d{1,2}(?=%)','match')));
  if usage >= maxusage
    fprintf(2,'%s: WARNING - PARTITION %s IS AT %d%% USAGE\n',...
      mfilename,partition,usage);
  end 
end
