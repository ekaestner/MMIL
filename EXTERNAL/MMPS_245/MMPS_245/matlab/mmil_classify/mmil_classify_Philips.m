function SeriesInfo = mmil_classify_Philips(SeriesInfo)
%function SeriesInfo = mmil_classify_Philips(SeriesInfo)
%
% Purpose: determine scan type for manufacturer Philips
%
% Required Input:
%   SeriesInfo: struct containing information about dicom series
%     initialized by mmil_serinfo
%
% Output:
%   SeriesInfo: struct containing series information and SeriesType
%
% Created:  05/30/11 by Don Hagler
% Last Mod: 12/21/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

for s=1:length(SeriesInfo)
  if SeriesInfo(s).ignore || ...
      isempty(regexp(upper(SeriesInfo(s).Manufacturer),'PHILIPS'))
    continue;
  end;
  if ~isempty(regexp(SeriesInfo(s).ScanningSequence,'GR')) &...
       (SeriesInfo(s).nimgs >= 124) & (SeriesInfo(s).SliceThickness <= 1.5) & ...
       (~isempty(regexp(upper(SeriesInfo(s).MRAcquisitionType),'3D')) |...
       ~isempty(regexp(upper(SeriesInfo(s).SequenceVariant),'MP')))
    SeriesInfo(s).SeriesType = 'MPR';
  elseif ~isempty(regexp(SeriesInfo(s).SeriesDescription,'SURVEY'))
    SeriesInfo(s).SeriesType = 'LOCALIZER';
  elseif ~isempty(regexp(SeriesInfo(s).ScanningSequence,'SE'))
    if ~isempty(regexp(SeriesInfo(s).SequenceName,'DWI')) |...
           ~isempty(regexp(SeriesInfo(s).SequenceName,'SEEPI'))
      if ~isempty(regexp(SeriesInfo(s).SeriesDescription,'DTI')) |...
             ~isempty(regexp(SeriesInfo(s).SeriesDescription,'Diff'))
         SeriesInfo(s).SeriesType = 'DTI';
         SeriesInfo(s).pepolar = 0;
      elseif SeriesInfo(s).NumberOfPhaseEncodingSteps >= 95 
         SeriesInfo(s).nb0 = 1;
         SeriesInfo(s).bval = 0;
         SeriesInfo(s).ndiffdirs = 0;
         if ~isempty(regexp(SeriesInfo(s).SeriesDescription,'B[0O]_A'))   
           SeriesInfo(s).SeriesType = 'DTI_rev';
           SeriesInfo(s).pepolar = 1;
         else
           SeriesInfo(s).SeriesType = 'DTI';
           SeriesInfo(s).pepolar = 0;
         end;
      else
        SeriesInfo(s).SeriesType = 'BOLD_SE'; 
        SeriesInfo(s).gradwarpinfo.unwarpflag = 0;        
        if ~isempty(regexp(SeriesInfo(s).SeriesDescription,'B[0O]_A'))   
          SeriesInfo(s).pepolar = 1;
        else
          SeriesInfo(s).pepolar = 0;
        end;
      end;               
    else 
      SeriesInfo(s).SeriesType = 'TSE2D'; % turbo spin echo
    end;
  elseif ~isempty(regexp(upper(SeriesInfo(s).SeriesDescription),'REST'))
    SeriesInfo(s).SeriesType = 'BOLD_bu';
    SeriesInfo(s).pepolar = 0;
  end;
end;

