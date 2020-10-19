function [ContainerOutInfo,ContainerOutDir] =...
    mmil_set_ContainerInfo(ContainerInInfo,ContainerInDir,ContainerType)
%function [ContainerOutInfo,ContainerOutDir] =...
%    mmil_set_ContainerInfo(ContainerInInfo,ContainerInDir,ContainerType)
%
% Required Input:
%   ContainerInInfo: ContainerInfo struct from input Container
%   ContainerInDir: name of input Container
%   ContainerType: type of output Container (e.g. MRIPROC, DTIPROC, etc.)
%
% Output:
%   ContainerOutInfo: ContainerInfo struct for output Container
%   ContainerOutDir: name of output Container
%
% Created:  09/11/12 by Don Hagler
% Last Mod: 09/11/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,3), return; end;

ContainerOutInfo = []; ContainerOutDir = [];

copy_tags = {'VisitID','StudyDate','StudyTime','StudyInstanceUID',...
  'MagneticFieldStrength','SeriesInfo','Manufacturer',...
  'ManufacturersModelName','MagneticFieldStrength'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set fields specific to ContainerOutInfo
ContainerOutInfo.SourceDir = ContainerInDir;
ContainerOutInfo.ContainerType = ContainerType;
ContainerOutInfo.ContainerUID = '1'; %% todo: system call to get unique container ID
ContainerOutInfo.ContainerCreationDate = datestr(now);

% copy fields from ContainerInInfo
for t=1:length(copy_tags)
  tag = copy_tags{t};
  ContainerOutInfo.(tag) = ContainerInInfo.(tag);
end;

% set output directory   
ContainerOutDir = sprintf('%s_%s_%s.%s_%s',ContainerType,...
  ContainerOutInfo.VisitID,ContainerOutInfo.StudyDate,...
  ContainerOutInfo.StudyTime,ContainerOutInfo.ContainerUID);

