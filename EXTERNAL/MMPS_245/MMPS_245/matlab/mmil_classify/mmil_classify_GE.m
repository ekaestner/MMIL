function SeriesInfo = mmil_classify_GE(SeriesInfo)
%function SeriesInfo = mmil_classify_GE(SeriesInfo)
%
% Purpose: determine scan type for manufacturer GE
%
% Required Input:
%   SeriesInfo: struct containing information about dicom series
%     initialized by mmil_serinfo
%
% Output:
%   SeriesInfo: struct containing series information and SeriesType
%
% Created:  05/30/11 by Don Hagler
% Last Mod: 02/04/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

for s=1:length(SeriesInfo)
  if SeriesInfo(s).ignore || ...
      isempty(regexp(upper(SeriesInfo(s).Manufacturer),'GE MEDICAL'))
    continue;
  end;
  invtime = SeriesInfo(s).InversionTime;
  magnitude_flag = SeriesInfo(s).magnitude_flag;
  % determine series type
  if magnitude_flag & SeriesInfo(s).FlipAngle>=10 &...
      ((~isempty(regexp(SeriesInfo(s).ScanOptions, 'FILTERED_GEMS')) &...
         ~isempty(regexp(SeriesInfo(s).SequenceName,'EFGRE3D')) & invtime>0) | ...
      (~isempty(regexp(SeriesInfo(s).SequenceName,'3DGRASS')) & invtime==0 &...
         SeriesInfo(s).RepetitionTime>10) | ...
      (~isempty(regexp(SeriesInfo(s).SeriesDescription,'FSPGR')) & invtime>0))
    SeriesInfo(s).SeriesType = 'FLASHhi';
  elseif magnitude_flag & invtime>0 &...
         ~isempty(regexp(SeriesInfo(s).SequenceName,'EFGRE3D'))
    SeriesInfo(s).SeriesType = 'MPR';
  elseif magnitude_flag & invtime==0 &...
         ~isempty(regexp(SeriesInfo(s).SequenceName,'3DGRASS')) &...
         SeriesInfo(s).FlipAngle<10 & SeriesInfo(s).RepetitionTime>10
    SeriesInfo(s).SeriesType = 'FLASHlo';
  elseif ~isempty(regexp(SeriesInfo(s).SequenceLabel,'xeta')) |...
         ~isempty(regexp(SeriesInfo(s).SequenceName,'spc3d1')) |...
         ~isempty(regexp(SeriesInfo(s).SequenceLabel,'Cube'))
    SeriesInfo(s).SeriesType = 'XetaT2';
  elseif magnitude_flag & ~isempty(regexp(SeriesInfo(s).SequenceName,'fl3d\d\d$'))
    SeriesInfo(s).SeriesType = 'T2STAR'; % DH - maybe not used
  elseif strcmp('epi2_pepolar',SeriesInfo(s).SequenceLabel)
    SeriesInfo(s).gradwarpinfo.unwarpflag = 0;
    SeriesInfo(s).SeriesType = 'DTI_ipp';
  elseif strcmp('epi2_pepolarFLEX',SeriesInfo(s).SequenceLabel)
    SeriesInfo(s).gradwarpinfo.unwarpflag = 0;
    if SeriesInfo(s).tensor_fnum > 0
      SeriesInfo(s).SeriesType = 'DTI_flex';
    else
      SeriesInfo(s).SeriesType = 'DTI_ipp';
    end;
  elseif magnitude_flag &...
         ~isempty(regexp(SeriesInfo(s).SeriesDescription,'DTI'))
    if ~isempty(regexp(upper(SeriesInfo(s).SequenceLabel),'MMILDTI')) % DH - RIL DTI
      SeriesInfo(s).gradwarpinfo.unwarpflag = 0;
    end;
    if ~isempty(regexp(upper(SeriesInfo(s).SeriesDescription),'REV'))
      SeriesInfo(s).SeriesType = 'DTI_rev';
    else
      SeriesInfo(s).SeriesType = 'DTI';
    end;
  elseif ~isempty(regexp(SeriesInfo(s).SequenceName,'fl3d1_ns$')) &...
         ismember(SeriesInfo(s).nimgs,[128 256]) &...
         SeriesInfo(s).Rows == 128 & SeriesInfo(s).Columns == 128 &...
         SeriesInfo(s).NumberOfPhaseEncodingSteps == 96
    SeriesInfo(s).SeriesType = 'AAS'; % DH - maybe not used
  elseif ~isempty(regexp(SeriesInfo(s).SeriesDescription,'Calibration Scan')) &...
         ~isempty(regexp(upper(SeriesInfo(s).MRAcquisitionType),'2D'))
    SeriesInfo(s).SeriesType = 'GEB1CAL'; % DH - calibration for PURE intensity normalization (done online on 1.5T at RIL)
  elseif ~isempty(regexp(upper(SeriesInfo(s).MRAcquisitionType),'2D')) & ...
         (~isempty(regexp(SeriesInfo(s).SequenceName,'FGRE')) |...
          ~isempty(regexp(SeriesInfo(s).SequenceName,'3.Plane')))
    SeriesInfo(s).SeriesType = 'LOCALIZER';
  elseif ~isempty(regexp(SeriesInfo(s).SequenceName,'EFGRE3D')) & SeriesInfo(s).RepetitionTime<10
    if strcmp(SeriesInfo(s).rfrxcoiltype,'UNKNOWN')
      fprintf('%s: WARNING: rfrxcoiltype = UNKNOWN\n',mfilename);
      SeriesInfo(s).SeriesType = 'UNKNOWN';
    elseif strcmp(SeriesInfo(s).rfrxcoiltype,'BODY')
      if magnitude_flag
        SeriesInfo(s).SeriesType = 'B1BC_mag'; % DH - old, only for ADNI, maybe VETSA1
      else
        SeriesInfo(s).SeriesType = 'B1BC_phi'; % DH - old, only for ADNI, maybe VETSA1
      end
    elseif strcmp(SeriesInfo(s).SeriesDescription,'T1 3D SAG Hi-Res')
      SeriesInfo(s).SeriesType = 'FLASHhi';
    else
      if magnitude_flag
        SeriesInfo(s).SeriesType = 'B1HC_mag'; % DH - old, only for ADNI, maybe VETSA1
      else
        SeriesInfo(s).SeriesType = 'B1HC_phi'; % DH - old, only for ADNI, maybe VETSA1
      end
    end
  elseif ~isempty(regexp(upper(SeriesInfo(s).MRAcquisitionType),'2D')) &...
         (~isempty(regexp(SeriesInfo(s).SequenceName,'FSE')) |...
          ~isempty(regexp(SeriesInfo(s).SeriesDescription,'FSE')))
    SeriesInfo(s).SeriesType = 'TSE2D'; % turbo spin echo
  elseif ~isempty(regexp(SeriesInfo(s).SequenceLabel,'epi_pepolar'))
    SeriesInfo(s).SeriesType = 'BOLD_pep';
    SeriesInfo(s).gradwarpinfo.unwarpflag = 0;
  elseif ~isempty(regexp(SeriesInfo(s).SequenceLabel,'epiRT_pepolarAPE'))
    SeriesInfo(s).SeriesType = 'BOLD_ape';
    SeriesInfo(s).gradwarpinfo.unwarpflag = 0;
  elseif ~isempty(regexp(SeriesInfo(s).ScanningSequence,'EP\\R'))
    if ~isempty(regexp(SeriesInfo(s).GE_SequenceLabel,'epiRT_pepolarIPP'))
      SeriesInfo(s).SeriesType = 'BOLD_ipp';
      SeriesInfo(s).gradwarpinfo.unwarpflag = 0;
    elseif ~isempty(regexp(SeriesInfo(s).GE_SequenceLabel,'epiRT_pepolar'))
      switch SeriesInfo(s).pepolar
        case 0
          SeriesInfo(s).SeriesType = 'BOLD_bu';
          % pepolar = 0: all frames are bottom-up ("forward")
        case 1
          SeriesInfo(s).SeriesType = 'BOLD_td';
          % pepolar = 1: all frames are top-down ("reverse")
        case {2,3}
          SeriesInfo(s).SeriesType = 'BOLD_ipp';
          % pepolar = 2: first frame is td, subsequent frames are bu
          % pepolar = 3: first frame is bu, subsequent frames are td
          SeriesInfo(s).gradwarpinfo.unwarpflag = 0; % grad warp is off for this sequence when pepolar = 2 or 3
        otherwise
          SeriesInfo(s).SeriesType = 'BOLD_bu';
      end;
    elseif ~isempty(regexp(SeriesInfo(s).SeriesDescription,'EPI_td_fmri')) |...
        ~isempty(regexp(SeriesInfo(s).SeriesDescription,'EPI td'))
      SeriesInfo(s).SeriesType = 'BOLD_td'; % DH - old data
    elseif ~isempty(regexp(SeriesInfo(s).SeriesDescription,'EPI_bu_fmri')) |...
        ~isempty(regexp(SeriesInfo(s).SeriesDescription,'EPI bu'))
      SeriesInfo(s).SeriesType = 'BOLD_bu'; % DH - old data
    else
     SeriesInfo(s).SeriesType = 'BOLD_bu';
    end;
  elseif ~isempty(regexp(SeriesInfo(s).ScanningSequence,'EP\\GR')) &...
        ~isempty(regexp(SeriesInfo(s).SequenceLabel,'epiRT'))
    if ~isempty(regexp(SeriesInfo(s).SeriesDescription,'fmri td')) |...
        ~isempty(regexp(SeriesInfo(s).SeriesDescription,'EPI td'))
      SeriesInfo(s).SeriesType = 'BOLD_td'; % DH - old data
    elseif ~isempty(regexp(SeriesInfo(s).SeriesDescription,'fmri bu')) |...
        ~isempty(regexp(SeriesInfo(s).SeriesDescription,'EPI bu'))
      SeriesInfo(s).SeriesType = 'BOLD_bu'; % DH - old data
    else
     SeriesInfo(s).SeriesType = 'BOLD_bu';
    end;
  elseif ~isempty(regexp(SeriesInfo(s).ScanningSequence,'EP\\G'))
    if ~isempty(regexp(SeriesInfo(s).SeriesDescription,'EPI_td_T1'))
      SeriesInfo(s).SeriesType = 'T1_EPI_td'; % DH - old data
    elseif ~isempty(regexp(SeriesInfo(s).SeriesDescription,'EPI_bu_T1'))
      SeriesInfo(s).SeriesType = 'T1_EPI_bu'; % DH - old data
    else
      SeriesInfo(s).SeriesType = 'EPI_misc';
    end;
  elseif ~isempty(regexp(SeriesInfo(s).SeriesDescription,'fm_TE1_NFS'))
    SeriesInfo(s).SeriesType = 'FMAP_TE1_NFS'; % field maps from cfMRI
  elseif ~isempty(regexp(upper(SeriesInfo(s).Manufacturer),'GE MEDICAL')) &...
      ~isempty(regexp(SeriesInfo(s).SeriesDescription,'fm_TE2_NFS'))
    SeriesInfo(s).SeriesType = 'FMAP_TE2_NFS'; % field maps from cfMRI
  elseif ~isempty(regexp(upper(SeriesInfo(s).MRAcquisitionType),'3D')) &...
         ~isempty(regexp(SeriesInfo(s).ScanningSequence,'GR')) &...
         SeriesInfo(s).nimgs >= 124 & SeriesInfo(s).SliceThickness <= 1.5
    % for 'SECONDARY' imgs from HAA_HUNT. AMD 12/09/10
    SeriesInfo(s).SeriesType = 'MPR';
  end;
end;

