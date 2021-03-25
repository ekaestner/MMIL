function [readoutFOV,phaseFOV,PhaseDir] = mmil_serinfo_FOV(serinfo)
%function [readoutFOV,phaseFOV,PhaseDir] = mmil_serinfo_FOV(serinfo)
%
% Required Input:
%   serinfo: series info struct created by mmil_classify_dicoms
%     must be a single struct, not an array
%
% Output:
%   readoutFOV: field of view for read-out direction
%   phaseFOV: field of view for phase-encode direction
%   PhaseDir: phase-encode direction (e.g. 'COL' or 'ROW')
%
% Created:  09/11/12 by Don Hagler
% Last Mod: 11/13/12 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
readoutFOV = []; phaseFOV = []; PhaseDir = [];

[nr,nc] = mmil_serinfo_dims(serinfo);
readoutFOV = serinfo.readoutFOV;
phaseFOV = serinfo.phaseFOV;
PhaseDir = upper(serinfo.InPlanePhaseEncodingDirection);

if isempty(readoutFOV)
  try
    if strcmp(PhaseDir,'COL')
      readoutFOV = serinfo.XPixelSpacing*nr;
      phaseFOV = serinfo.YPixelSpacing*nc;
    elseif strcmp(PhaseDir,'ROW')
      readoutFOV = serinfo.YPixelSpacing*nc;
      phaseFOV = serinfo.XPixelSpacing*nr;
    else
      fprintf('%s: WARNING: PhaseDir = %s\n',mfilename,PhaseDir);
      readoutFOV = serinfo.YPixelSpacing*nc;
      phaseFOV = serinfo.XPixelSpacing*nr;
    end;
  catch
    fprintf('%s: WARNING: failed to get FOV:\n%s\n',mfilename,lasterr);
  end;
end;

