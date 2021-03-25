function SeriesInfo = mmil_qmat_DTI(SeriesInfo)
%function SeriesInfo = mmil_qmat_DTI(SeriesInfo)
%
% Purpose: modify SeriesInfo by setting qmat, nb0, ndiffdirs, and bval
%   for series with 'DTI' in the SeriesType
%
% Required Input:
%   SeriesInfo: struct array containing information about
%     multiple dicom image series including SeriesType
%
% Output:
%   SeriesInfo: struct array with slice timing information for BOLD series
%
% Created:  11/14/12 by Don Hagler
% Last Mod: 11/16/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

for s=1:length(SeriesInfo)
  if SeriesInfo(s).ignore, continue; end;
  if strncmp(SeriesInfo(s).SeriesType,'DTI',3)
    if mmil_Philips(SeriesInfo(s))
      if SeriesInfo(s).nb0~=1
        [SeriesInfo(s).diffdirs,SeriesInfo(s).nb0,SeriesInfo(s).ndiffdirs,...
         SeriesInfo(s).bval] = mmil_qmat_Philips(SeriesInfo(s));
      end;
    elseif mmil_GE(SeriesInfo(s))
      if isfield(SeriesInfo(s).info,'Private_0019_10bb') &&...
         isfield(SeriesInfo(s).info,'Private_0019_10bc') &&...
         isfield(SeriesInfo(s).info,'Private_0019_10bd') &&...
         isfield(SeriesInfo(s).info,'Private_0043_1039')
        [SeriesInfo(s).diffdirs,SeriesInfo(s).nb0,ndd,SeriesInfo(s).bval] =...
           mmil_qmat_GE(SeriesInfo(s));  
        if ndd && ndd ~= SeriesInfo(s).ndiffdirs
          fprintf('%s: WARNING: ndd (%d) and ndiffdirs (%d) do not match for Series(%d)\n',...
            mfilename,ndd,SeriesInfo(s).ndiffdirs,s);
          SeriesInfo(s).ndiffdirs = ndd;
        end; 
      end;
    elseif mmil_Siemens(SeriesInfo(s))
      if isfield(SeriesInfo(s).info, 'Private_0019_100c')
        if ~isempty(regexp(SeriesInfo(s).ImageType, 'MOSAIC'))
          [SeriesInfo(s).diffdirs,SeriesInfo(s).nb0,ndd,bval] = ...
            mmil_qmat_Siemens(SeriesInfo(s));
          if ndd ~= SeriesInfo(s).ndiffdirs
            fprintf('%s: WARNING: ndd (%d) and ndiffdirs (%d) do not match\n',...
              mfilename,ndd,SeriesInfo(s).ndiffdirs);
            SeriesInfo(s).ndiffdirs = ndd;
          end;
          if SeriesInfo(s).bval~=bval
            fprintf('%s: WARNING: bval was %d, replaced with %d\n',...
              mfilename,SeriesInfo(s).bval,bval);
            SeriesInfo(s).bval = bval;
          end;
        else %% todo: what if actual DTI scan is non-mosaic?
          SeriesInfo(s).nb0 = 1;
          SeriesInfo(s).ndiffdirs = 0;
        end; 
      end;
      if SeriesInfo(s).nb0 == 0
        fprintf('%s: WARNING: nb0 = 0, replacing with nb0 = 1\n',mfilename);
        SeriesInfo(s).nb0 = 1;
      end
    end;
  end;
end;

