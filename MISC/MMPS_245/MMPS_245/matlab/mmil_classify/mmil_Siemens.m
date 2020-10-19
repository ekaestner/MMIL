function res = mmil_Siemens(mri_info)
%function res = mmil_Siemens(mri_info)
%
% Required Input:
%   mri_info: struct or struct array containing field 'Manufacturer'
%     may be ContainerInfo, SeriesInfo, serinfo, etc.
%
% Output:
%   res: 1 if mri_info.Manufacturer contains SIEMENS
%        0 if not
%
% Created:  09/11/12 by Don Hagler
% Last Mod: 12/03/12 by Don Hagler
%

res = 0;
if ~isfield(mri_info,'Manufacturer')
  error('mri_info is missing Manufacturer field');
else
  Manufacturer = upper({mri_info.Manufacturer});
  if any(~cellfun(@isempty,regexp(Manufacturer,'SIEMENS')))
    res = 1;
  end;
end;

