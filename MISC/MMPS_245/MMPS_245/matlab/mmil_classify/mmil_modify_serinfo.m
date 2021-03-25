function SeriesInfo = mmil_modify_serinfo(SeriesInfo)
%function SeriesInfo = mmil_modify_serinfo(SeriesInfo)
%
% Purpose: modify SeriesInfo by setting fields based on other information
%   inlcluding manufacturer specific private tags
%
% Required Input:
%   SeriesInfo: struct array containing information about
%     multiple dicom image series initialized by mmil_init_serinfo
%
% Output:
%   SeriesInfo: struct array with modified information
%
% Created:  11/13/12 by Don Hagler
% Last Mod: 02/04/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

for s=1:length(SeriesInfo)
  fname = SeriesInfo(s).FileNames{1};
  if SeriesInfo(s).ignore, continue; end;
  if strcmp(SeriesInfo(s).ManufacturersModelName,'UNDEFINED') &&...
     isfield(SeriesInfo(s).info,'ManufacturerModelName')
    SeriesInfo(s).ManufacturersModelName = ...
      SeriesInfo(s).info.ManufacturerModelName;
  end;
  if length(SeriesInfo(s).PixelSpacing)==2
    SeriesInfo(s).XPixelSpacing = SeriesInfo(s).PixelSpacing(1);
    SeriesInfo(s).YPixelSpacing = SeriesInfo(s).PixelSpacing(2);
  end;
  [SeriesInfo(s).gradwarpinfo,errmsg] = ctx_get_gradwarpinfo(SeriesInfo(s).info);
  if ~isempty(errmsg)
    fprintf('%s: WARNING: failed to get gradwarpinfo for %s:\n%s\n',...
      mfilename,fname,errmsg);
  end
  if ~isempty(SeriesInfo(s).gradwarpinfo) && ...
     ~isfield(SeriesInfo(s).gradwarpinfo,'gwtype')
    fprintf('%s: WARNING: missing gwtype for %s\n',mfilename,fname);
  end
  if ~isempty(SeriesInfo(s).gradwarpinfo) &&...
      isfield(SeriesInfo(s).gradwarpinfo,'unwarpflag') && ...
      SeriesInfo(s).gradwarpinfo.unwarpflag
    % determine whether grad warp was really done online
    n = 1+floor(SeriesInfo(s).nimgs/2);
    unwarpflag = mmil_detect_gradwarp(SeriesInfo(s).FileNames{n});
    if ~unwarpflag, SeriesInfo(s).gradwarpinfo.unwarpflag = 0; end;
  end;
end;

% get additional information from manufacturer-specific private tags
if mmil_Philips(SeriesInfo)
  SeriesInfo = mmil_modify_serinfo_Philips(SeriesInfo);
elseif mmil_GE(SeriesInfo)
  SeriesInfo = mmil_modify_serinfo_GE(SeriesInfo);
elseif mmil_Siemens(SeriesInfo)
  SeriesInfo = mmil_modify_serinfo_Siemens(SeriesInfo);
end;

SeriesInfo = mmil_modify_serinfo_rfrxcoiltype(SeriesInfo);

