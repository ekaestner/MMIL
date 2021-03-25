function [ContainerOutDir,errcode] = MMIL_Convert_MRIRAW_to_MRIPROC(...
  ContainerRootInDir,ContainerInDir,ContainerRootOutDir,forceflag)
%function [ContainerOutDir,errcode] = MMIL_Convert_MRIRAW_to_MRIPROC(...
%  ContainerRootInDir,ContainerInDir,ContainerRootOutDir,forceflag)
%
% Required Input:
%   ContainerRootInDir: path of root directory for MRIRAW containers
%   ContainerInDir: input MRIRAW container
%   ContainerRootOutDir: path of root dirrectory for MRIPROC containers
%
% Optional Input:
%   forceflag: [0|1] overwrite existing output
%     {default = 0}
%
% Created:  12/16/08 by Don Hagler
% Prev Mod: 04/22/16 by Sean Hatton
% Last Mod: 04/22/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ContainerOutDir = []; errcode = 0;

if ~mmil_check_nargs(nargin,3), return; end;
if ~exist('forceflag','var'), forceflag = false; end;

% must have one of these
required_SeriesTypes = ...
  {'MPR','FLASHhi','FLASHlo','XetaT2','MEDIChi','MEDIClo'};

% may have any of these
valid_SeriesTypes = union(required_SeriesTypes,...
  {'B1HC_mag','B1HC_phi','B1BC_mag','B1BC_phi','GEB1CAL','TSE2D'});

% fields to be initialized in ScanInfo
output_ScanTypes = {'MPR','FLASHhi','FLASHlo','XetaT2',...
  'MEDIChi','MEDIClo','B1HC','B1BC','GEB1CAL','TSE2D'};

fprintf('%s(''%s'',''%s'',''%s'',%d)\n',...
  mfilename,ContainerRootInDir,ContainerInDir,ContainerRootOutDir,forceflag);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get input ContainerInfo
ContainerInPath = sprintf('%s/%s',ContainerRootInDir,ContainerInDir);
[ContainerInInfo,errcode] = MMIL_Load_ContainerInfo(ContainerInPath);
if errcode, return; end;

% set output ContainerInfo, create ContainerOutPath
[ContainerOutInfo,ContainerOutDir] =...
  mmil_set_ContainerInfo(ContainerInInfo,ContainerInDir,'MRIPROC');
ContainerOutPath = sprintf('%s/%s',ContainerRootOutDir,ContainerOutDir);

% check if there are any of the required series in this study, return if not
SeriesInfo = ContainerInInfo.SeriesInfo;
if isempty(intersect({SeriesInfo.SeriesType},required_SeriesTypes))
  fprintf('%s: WARNING: no structural MRI scans to convert\n',mfilename);
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
    case {'MPR','FLASHhi','FLASHlo','XetaT2','TSE2D'}
      ScanInfo.(serinfo.SeriesType) = convert_structural(serinfo,s,...
        ContainerOutPath,ScanInfo.(serinfo.SeriesType),forceflag);
    case {'MEDIChi','MEDIClo'}
      ScanInfo.(serinfo.SeriesType) = convert_MEDIC(serinfo,s,...
        ContainerOutPath,ScanInfo.(serinfo.SeriesType),forceflag);
    case {'B1HC_mag','B1HC_phi'}
      ScanInfo.B1HC = ...
        convert_B1(serinfo,s,ContainerOutPath,ScanInfo.B1HC,forceflag);
    case {'B1BC_mag','B1BC_phi'}
      ScanInfo.B1BC = ...
        convert_B1(serinfo,s,ContainerOutPath,ScanInfo.B1BC,forceflag);
    case 'GEB1CAL'
      ScanInfo.GEB1CAL = ...
        convert_GEB1CAL(serinfo,s,ContainerOutPath,ScanInfo.GEB1CAL,forceflag);
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

function all_scaninfo = convert_structural(serinfo,i_series,outdir,...
                                    all_scaninfo,forceflag)
  if mmil_check_nimgs_mismatch(serinfo,i_series), return; end;
  i_scan = length(all_scaninfo) + 1;
  scaninfo = mmil_convert_structural(serinfo,i_series,i_scan,outdir,forceflag);
  if isempty(all_scaninfo)
    all_scaninfo = scaninfo;
  else
    all_scaninfo(end+1) = scaninfo;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function all_scaninfo = convert_MEDIC(serinfo,i_series,outdir,...
                                    all_scaninfo,forceflag)
  i_scan = length(all_scaninfo) + 1;
  scaninfo = mmil_convert_MEDIC(serinfo,i_series,i_scan,outdir,forceflag);
  if isempty(all_scaninfo)
    all_scaninfo = scaninfo;
  else
    all_scaninfo(end+1) = scaninfo;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function all_scaninfo = convert_B1(serinfo,i_series,outdir,...
                                    all_scaninfo,forceflag)
  if mmil_check_nimgs_mismatch(serinfo,i_series), return; end;
  i_scan = length(all_scaninfo) + 1;
  scaninfo = mmil_convert_B1(serinfo,i_series,i_scan,outdir,forceflag);
  if isempty(all_scaninfo)
    all_scaninfo = scaninfo;
  else
    all_scaninfo(end+1) = scaninfo;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function all_scaninfo = convert_GEB1CAL(serinfo,i_series,outdir,...
                                    all_scaninfo,forceflag)
  if mmil_check_nimgs_mismatch(serinfo,i_series), return; end;
  i_scan = length(all_scaninfo) + 1;
  scaninfo = mmil_convert_GEB1CAL(serinfo,i_series,i_scan,outdir,forceflag);
  if isempty(all_scaninfo)
    all_scaninfo = scaninfo;
  else
    all_scaninfo(end+1) = scaninfo;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
