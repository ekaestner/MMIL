function [gradwarpinfo,errmsg] = ctx_get_gradwarpinfo(dcminfo)
%function [gradwarpinfo,errmsg] = ctx_get_gradwarpinfo(dcminfo)
%
% Early Mod: 05/19/10 by Don Hagler
% Last Mod:  01/21/17 by Don Hagler
%

errmsg='';

Manufacturer = dcminfo(1).Manufacturer;
if isfield(dcminfo(1),'ManufacturersModelName')
  ManufacturersModelName = dcminfo(1).ManufacturersModelName;
else
  ManufacturersModelName = dcminfo(1).ManufacturerModelName;
end
switch lower(Manufacturer)
  case 'siemens'
    gradwarpinfo.unwarpflag = 0; % full 3D unwarping
    gradwarpinfo.isoctrflag = 0; % assume not using isocenter scanning (NB: should check DICOM header info on new systems)
    switch lower(ManufacturersModelName)
      case {'sonata','trio','sonatavision'}
        gradwarpinfo.gwtype = 0;
      case {'allegra'}
        gradwarpinfo.gwtype = 1;
      case {'avanto','triotim'}
        gradwarpinfo.gwtype = 4;
      case {'espree','axxess'}
        gradwarpinfo.gwtype = 5;
      case {'symphony','symphonyvision'} % Quantum
        gradwarpinfo.gwtype = 6;
      case {'skyra'} %% todo: check this
        gradwarpinfo.gwtype = 11;
      case {'connectome'} %% todo: check this
        gradwarpinfo.gwtype = 12;
      case {'prisma','prisma_fit'}
        gradwarpinfo.gwtype = 13;
      otherwise
        errmsg=sprintf('%s: Unknown gradient model %s %s\n',mfilename,Manufacturer,ManufacturersModelName);
        return;
    end
  case 'ge medical systems'
    gradwarpinfo.unwarpflag = 1; % through-plane only (should check NOGRAD CV)
    gradwarpinfo.isoctrflag = 1; % assume using isocenter scanning (should check DICOM tag?)
    switch deblank(lower(ManufacturersModelName))
      case {'discovery mr450','discovery mr750'}
        gradwarpinfo.gwtype = 9;
        return;
      case {'discovery mr750w'} %% todo: check this
        gradwarpinfo.gwtype = 10;
        return;
    end
    if isfield(dcminfo(1),'Private_0043_106f')
      tmp = dcminfo(1).Private_0043_106f;
      if length(tmp)>=2 & tmp(1)>=48 & tmp(1)<=57 & tmp(2)==92 % Does it look like ASCII nums delimited by \?
        tmp = tmp(1:2:end)-48;
      end
      if length(tmp) >=4
        key = tmp(4);
      else
        key = 0;
      end
      if key == 1
        gradwarpinfo.gwtype = 7; % Whole
      elseif key == 2;
        gradwarpinfo.gwtype = 8; % Zoom
      elseif key == 0;
%        gradwarpinfo.gwtype = 2; % BRM mode??? NB: Needs to be checked!
        gradwarpinfo.ambiguousgwtype = 1; % Ambiguous
        errmsg=sprintf('%s: Ambiguous gradient info for %s %s\n',mfilename,Manufacturer,ManufacturersModelName);
        return;
      else
        gradwarpinfo.ambiguousgwtype = 1; % Ambiguous
        errmsg=sprintf('%s: TwinSpeed mode = %d unknown for %s %s\n',mfilename,key,Manufacturer,ManufacturersModelName);
        return;
      end
      return;
    end
    switch lower(ManufacturersModelName)
      case {'brm'} % Check on actual name
        gradwarpinfo.gwtype = 2;
      case {'crm'} % Check on actual name
        gradwarpinfo.gwtype = 3;
      case {'signa excite','signa hdx'}
        gradwarpinfo.ambiguousgwtype = 1; % Ambiguous
        errmsg=sprintf('%s: Missing DICOM tag Private_0043_106f for %s %s\n',mfilename,Manufacturer,ManufacturersModelName);
        return;
      case 'genesis_signa'
%        gradwarpinfo.gwtype = 2; % BRM mode??? NB: Needs to be checked!
%        GET INFO FROM TEXTFILE FOR SYSTEM           XXX
        gradwarpinfo.ambiguousgwtype = 1; % Ambiguous
        errmsg=sprintf('%s: Ambiguous gradient info for %s %s\n',mfilename,Manufacturer,ManufacturersModelName);
        return;
      otherwise
        gradwarpinfo.ambiguousgwtype = 1; % Ambiguous
        errmsg=sprintf('%s: Unknown gradient model %s %s\n',mfilename,Manufacturer,ManufacturersModelName);
        return;
    end
  case 'philips medical systems'
    gradwarpinfo = []; % Default to no gradwarp for Philips scanners
  otherwise
    gradwarpinfo = []; % Default to no gradwarp for unknown scanner mfgrs
  end

