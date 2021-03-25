function scaninfo = mmil_convert_BOLD(serinfo,i_series,i_scan,outdir,forceflag)
%function scaninfo = mmil_convert_BOLD(serinfo,i_series,i_scan,outdir,forceflag)
%
% Required Input:
%   serinfo: series info struct created by mmil_classify_dicoms
%     must be a single struct, not an array
%   i_series: series index (within study)
%   i_scan: scan index (used in output name)
%
% Optional Input:
%   outdir: output directory
%     {default = pwd}
%   forceflag: [0|1] overwrite existing output
%     {default = 0}
%
% Output:
%   scaninfo: struct containing selected information for scan
%
% Created:  02/03/16 by Don Hagler
% Last Mod: 02/03/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scaninfo = [];
if ~mmil_check_nargs(nargin,3), return; end;
if ~exist('outdir','var') || isempty(outdir), outdir = pwd; end;
if ~exist('forceflag','var') || isempty(forceflag), forceflag = 0; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

