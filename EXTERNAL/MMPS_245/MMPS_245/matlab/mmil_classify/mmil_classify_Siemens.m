function SeriesInfo = mmil_classify_Siemens(SeriesInfo)
%function SeriesInfo = mmil_classify_Siemens(SeriesInfo)
%
% Purpose: determine scan type for manufacturer Siemens
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
      isempty(regexp(upper(SeriesInfo(s).Manufacturer),'SIEMENS'))
    continue;
  end;
  SeriesInfo(s).nb0 = 1;  % Only used for DTI
  magnitude_flag = SeriesInfo(s).magnitude_flag;
  if (magnitude_flag & ~isempty(regexp(SeriesInfo(s).SequenceName,'tfl3d1'))) | ...
     (~isempty(regexp(SeriesInfo(s).SeriesDescription,'mpr3D'))) | ... % for TOP (excessive anonymization) 06/25/08 DH
     (~isempty(regexp(SeriesInfo(s).SequenceName,'tfl3d1')) &...
      ~isempty(findstr('ORIGINAL\PRIMARY\M\DIS2D',SeriesInfo(s).ImageType)))
      SeriesInfo(s).SeriesType = 'MPR';
  elseif ~isempty(regexp(SeriesInfo(s).SequenceName,'fl3d\d_ns')) &...
         SeriesInfo(s).FlipAngle>=10 & SeriesInfo(s).RepetitionTime>10 &...
         magnitude_flag
    SeriesInfo(s).SeriesType = 'MEDIChi'; % DH - old, not used since MGH (VETSA1)
  elseif ~isempty(regexp(SeriesInfo(s).SequenceName,'fl3d\d_ns')) &...
         SeriesInfo(s).FlipAngle<10 & SeriesInfo(s).RepetitionTime>10 &...
         magnitude_flag
    SeriesInfo(s).SeriesType = 'MEDIClo'; % DH - old, not used since MGH (VETSA1)
  elseif ~isempty(regexp(SeriesInfo(s).SequenceName,'fl3d\d')) &...
         SeriesInfo(s).FlipAngle<10 & SeriesInfo(s).RepetitionTime>10 &...
         magnitude_flag
    SeriesInfo(s).SeriesType = 'FLASHlo'; % DH - not commonly used (only for Thomas Thesen)
  elseif ~isempty(regexp(SeriesInfo(s).SequenceName,'fl3d\d')) &...
         SeriesInfo(s).FlipAngle>=10 & SeriesInfo(s).RepetitionTime>10 &...
         magnitude_flag
    SeriesInfo(s).SeriesType = 'FLASHhi'; % DH - not commonly used (only for Thomas Thesen)
  elseif ~isempty(regexp(SeriesInfo(s).SequenceName,'spc3d1')) 
    SeriesInfo(s).SeriesType = 'XetaT2';
  elseif ~isempty(regexp(SeriesInfo(s).SequenceName,'epse2d1_96'))
    SeriesInfo(s).nb0 = 1;
    if SeriesInfo(s).pe_extra_val > 10
      SeriesInfo(s).SeriesType = 'DTI_rev';
      SeriesInfo(s).pepolar = 1;
    else
      SeriesInfo(s).SeriesType = 'DTI';
      SeriesInfo(s).pepolar = 0;
    end;
  elseif ~isempty(regexp(SeriesInfo(s).SequenceName,'epse2d1_64'))
    SeriesInfo(s).nb0 = 1;
    if SeriesInfo(s).nimgs == 2
      % if 2 images, then it is a mosaic field map,
      %   with each image containing a different phase encoding direction
      SeriesInfo(s).SeriesType = 'BOLD_pep';
      if SeriesInfo(s).pe_extra_val > 4
        SeriesInfo(s).pepolar = 2;
      else
        SeriesInfo(s).pepolar = 3;
      end;
    else
      SeriesInfo(s).SeriesType = 'BOLD_SE';
      if SeriesInfo(s).pe_rev_flag | (SeriesInfo(s).pe_extra_val > 4)
        SeriesInfo(s).pepolar = 1;
      else
        SeriesInfo(s).pepolar = 0;
      end;
      SeriesInfo(s).gradwarpinfo.unwarpflag = 0;
    end;
  elseif ~isempty(regexp(SeriesInfo(s).SequenceName,'ep_b')) &...
         (~isempty(regexp(upper(SeriesInfo(s).SeriesDescription),'P-A')) | ...
           SeriesInfo(s).pe_rev_flag |...
           (~isempty(SeriesInfo(s).InPlaneRot) && ...
              abs(SeriesInfo(s).InPlaneRot)>0.2))
    if ~isempty(regexp(upper(SeriesInfo(s).SeriesDescription),'BOLD')) |...
        SeriesInfo(s).NumberOfPhaseEncodingSteps <= 64
      SeriesInfo(s).SeriesType = 'BOLD_SE'; 
      SeriesInfo(s).gradwarpinfo.unwarpflag = 0;        
      SeriesInfo(s).pepolar = 1;
    else
      SeriesInfo(s).SeriesType = 'DTI_rev';
      SeriesInfo(s).pepolar = 1;
    end;
  elseif ~isempty(regexp(SeriesInfo(s).SequenceName,'ep_b'))  &...
         (~isempty(regexp(upper(SeriesInfo(s).SeriesDescription),'A-P')) | ...
           ~isempty(regexp(upper(SeriesInfo(s).SeriesDescription),'AP')) | ...
           SeriesInfo(s).pe_for_flag) & isempty(regexp(upper(SeriesInfo(s).SeriesDescription),'P-A'))
    if ~isempty(regexp(upper(SeriesInfo(s).SeriesDescription),'BOLD')) |...
       SeriesInfo(s).NumberOfPhaseEncodingSteps <= 64
      SeriesInfo(s).SeriesType = 'BOLD_SE'; 
      SeriesInfo(s).gradwarpinfo.unwarpflag = 0;        
      SeriesInfo(s).pepolar = 0;
    else
      SeriesInfo(s).SeriesType = 'DTI';
      SeriesInfo(s).pepolar = 0;
    end;
  elseif ~isempty(regexp(SeriesInfo(s).SequenceName,'ep_b')) & magnitude_flag &...
         (~isempty(regexp(SeriesInfo(s).SeriesDescription,'rorevro')) |...
           ~isempty(regexp(SeriesInfo(s).SeriesDescription,'DTI')))
    SeriesInfo(s).SeriesType = 'DTI';
    SeriesInfo(s).pepolar = 3;
  elseif (~isempty(regexp(upper(SeriesInfo(s).MRAcquisitionType),'2D')) &...
         ~isempty(regexp(SeriesInfo(s).SequenceName,'fl2d1')) &...
         SeriesInfo(s).nimgs <= 6)
    SeriesInfo(s).SeriesType = 'LOCALIZER';
  elseif (~isempty(regexp(SeriesInfo(s).SequenceName,'fl3d1_ns$')) &...
         SeriesInfo(s).Rows == 128 & SeriesInfo(s).Columns == 128 &...
         SeriesInfo(s).NumberOfPhaseEncodingSteps == 128)
    if strcmp(SeriesInfo(s).rfrxcoiltype,'UNKNOWN')
      fprintf('%s: WARNING: rfrxcoiltype = UNKNOWN\n',mfilename);
      SeriesInfo(s).SeriesType = 'UNKNOWN';
    elseif strcmp(SeriesInfo(s).rfrxcoiltype,'BODY')
      if magnitude_flag
        SeriesInfo(s).SeriesType = 'B1BC_mag';
      else
        SeriesInfo(s).SeriesType = 'B1BC_phi';
      end
    else
      if magnitude_flag
        SeriesInfo(s).SeriesType = 'B1HC_mag';
      else
        SeriesInfo(s).SeriesType = 'B1HC_phi';
      end
    end
  elseif ~isempty(regexp(SeriesInfo(s).SequenceName,'tse2d'))
    SeriesInfo(s).SeriesType = 'TSE2D'; % turbo spin echo
  elseif ~isempty(regexp(SeriesInfo(s).SeriesDescription,...
                         'gre_mgh_me_B0_field_map'))
    SeriesInfo(s).SeriesType = 'B0FM';
  elseif ~isempty(regexp(SeriesInfo(s).SequenceName,'fm2d2r'))
    %SeriesInfo(s).SeriesType = 'FMAP'; 
    SeriesInfo(s).SeriesType = 'JUNK';
%  elseif ~isempty(regexp(SeriesInfo(s).SequenceName,'epfid2d1_64'))
  elseif ~isempty(regexp(SeriesInfo(s).SequenceName,'epfid2d1'))
    SeriesInfo(s).gradwarpinfo.unwarpflag = 0;
    if SeriesInfo(s).pe_rev_flag
      SeriesInfo(s).SeriesType = 'BOLD_td';
      SeriesInfo(s).pepolar = 1;
    else
      SeriesInfo(s).SeriesType = 'BOLD_bu';
      SeriesInfo(s).pepolar = 0;
    end;
  elseif ~isempty(regexp(SeriesInfo(s).SequenceName,'tfi2d1'))
    SeriesInfo(s).SeriesType = 'JUNK';
  end;
end;

