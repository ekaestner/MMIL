function [ContainerOutDir,errcode] = MMIL_Convert_MRIRAW_to_BOLDPROC(...
  ContainerRootInDir,ContainerInDir,ContainerRootOutDir,forceflag)
%function [ContainerOutDir,errcode] = MMIL_Convert_MRIRAW_to_BOLDPROC(...
%  ContainerRootInDir,ContainerInDir,ContainerRootOutDir,forceflag)
%
% Required Input:
%   ContainerRootInDir: path of root directory for MRIRAW containers
%   ContainerInDir: input MRIRAW container
%   ContainerRootOutDir: path of root dirrectory for BOLDPROC containers
%
% Optional Input:
%   forceflag: [0|1] overwrite existing output
%     {default = 0}
%
% Created:  09/10/12 by Don Hagler
% Last Mod: 03/20/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ContainerOutDir = []; errcode = 0;

if ~mmil_check_nargs(nargin,3), return; end;
if ~exist('forceflag','var'), forceflag = false; end;

% must have one of these
required_SeriesTypes = ...
  {'BOLD','BOLD_td','BOLD_bu','BOLD_ipp','BOLD_pep','BOLD_SE','BOLD_ape'};

% may have any of these
valid_SeriesTypes = union(required_SeriesTypes,...
  {'FMAP_TE1_NFS','FMAP_TE2_NFS','FMAP'});

% fields to be initialized in ScanInfo
output_ScanTypes = {'BOLD','FMAP'};

fprintf('%s(''%s'',''%s'',''%s'',%d)\n',...
  mfilename,ContainerRootInDir,ContainerInDir,ContainerRootOutDir,forceflag);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get input ContainerInfo
ContainerInPath = sprintf('%s/%s',ContainerRootInDir,ContainerInDir);
[ContainerInInfo,errcode] = MMIL_Load_ContainerInfo(ContainerInPath);
if errcode, return; end;

% set output ContainerInfo, create ContainerOutPath
[ContainerOutInfo,ContainerOutDir] =...
  mmil_set_ContainerInfo(ContainerInInfo,ContainerInDir,'BOLDPROC');
ContainerOutPath = sprintf('%s/%s',ContainerRootOutDir,ContainerOutDir);

% check if there are any of the required series in this study, return if not
SeriesInfo = ContainerInInfo.SeriesInfo;
if isempty(intersect({SeriesInfo.SeriesType},required_SeriesTypes))
  fprintf('%s: WARNING: no BOLD scans to convert\n',mfilename);
  errcode = 1;
  return;
end;
mmil_mkdir(ContainerOutPath);

ScanInfo = mmil_init_ScanInfo(output_ScanTypes);
for s=1:length(SeriesInfo)
  % check series type, get info, fnames
  [serinfo,valid_flag,missing_flag] = ...
    mmil_check_series(SeriesInfo,s,ContainerInPath,valid_SeriesTypes);
  if ~valid_flag || missing_flag, continue; end;
  switch serinfo.SeriesType
    case {'BOLD','BOLD_td','BOLD_bu','BOLD_ipp','BOLD_pep','BOLD_SE','BOLD_ape'}
      ScanInfo.BOLD = ...
        convert_BOLD(serinfo,s,ContainerOutPath,ScanInfo.BOLD,forceflag);
    case {'FMAP_TE1_NFS','FMAP_TE2_NFS','FMAP'}
      ScanInfo.FMAP = ...
        convert_FMAP(serinfo,s,ContainerOutPath,ScanInfo.FMAP,forceflag);
  end
end

% save ContainerOutInfo to ContainerInfo.mat file
errcode = MMIL_Save_ContainerInfo_ScanInfo(ContainerOutPath,...
                                           ContainerOutInfo,ScanInfo);
if errcode, return; end;

% save SeriesIndex to ScanInfo.csv file
mmil_write_ScanInfo(ScanInfo,ContainerOutPath);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function all_scaninfo = convert_BOLD(serinfo,i_series,outdir,all_scaninfo,forceflag)
  i_scan = length(all_scaninfo) + 1;
  if mmil_Philips(serinfo)
    scaninfo = mmil_convert_BOLD_Philips(serinfo,i_series,i_scan,outdir,forceflag);
  elseif mmil_Siemens(serinfo)
    scaninfo = mmil_convert_BOLD_Siemens(serinfo,i_series,i_scan,outdir,forceflag);
  elseif mmil_GE(serinfo)
    scaninfo = mmil_convert_BOLD_GE(serinfo,i_series,i_scan,outdir,forceflag);
  else
    fprintf('%s: %s is unsupported manufacturer for BOLD type\n',...
      mfilename,serinfo.Manufacturer);
    return;
  end;
  for i=1:length(scaninfo)
    if isempty(all_scaninfo)
      all_scaninfo = scaninfo(i);
    else
      all_scaninfo(end+1) = scaninfo(i);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function all_scaninfo = convert_FMAP(serinfo,i_series,outdir,all_scaninfo,forceflag)
  i_scan = length(all_scaninfo) + 1;
  if mmil_Siemens(serinfo) || mmil_GE(serinfo)
    scaninfo = mmil_convert_FMAP(serinfo,i_series,i_scan,outdir,forceflag);
  else
    fprintf('%s: %s is unsupported manufacturer for FMAP type\n',...
      mfilename,serinfo.Manufacturer);
    return;
  end;
  for i=1:length(scaninfo)
    if isempty(all_scaninfo)
      all_scaninfo = scaninfo(i);
    else
      all_scaninfo(end+1) = scaninfo(i);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
