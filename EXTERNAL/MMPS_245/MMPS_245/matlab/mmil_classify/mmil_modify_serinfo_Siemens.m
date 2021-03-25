function SeriesInfo = mmil_classify_Siemens(SeriesInfo)
%function SeriesInfo = mmil_classify_Siemens(SeriesInfo)
%
% Purpose: gather information from dicom header for manufacturer Siemens
%
% Required Input:
%   SeriesInfo: struct containing information about dicom
%     initialized by mmil_classify_dicom
%
% Output:
%   SeriesInfo: struct containing additional information
%     extracted from Siemens-specific private tags
%
% Created:  11/13/12 by Don Hagler
% Last Mod: 04/01/15 by Josh Kuperman
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

for s=1:length(SeriesInfo)
  if SeriesInfo(s).ignore, continue; end;
  if isfield(SeriesInfo(s).info,'Private_0029_2000')
    tmp = char(SeriesInfo(s).info.Private_0029_2000);
    tmp = reshape(tmp,[1,prod(size(tmp))]);
    SeriesInfo(s).scalefact = str2double(tmp);
  end;
  if isfield(SeriesInfo(s).info,'Private_0029_1020')
    [mdhStruct,txt] = mmil_read_mdhStruct(SeriesInfo(s).info);
    if isfield(mdhStruct,...
              'asCoilSelectMeasl0rdasListl0rdsCoilElementIDdtCoilID')
      SeriesInfo(s).ReceiveCoilName = ...
        regexprep(mdhStruct.asCoilSelectMeasl0rdasListl0rdsCoilElementIDdtCoilID,...
          '"','');;
    elseif isfield(mdhStruct,...
                  'asCOIL_SELECT_MEASl0rdasListl0rdsCoilElementIDdtCoilID')
      SeriesInfo(s).ReceiveCoilName = ...
        regexprep(mdhStruct.asCOIL_SELECT_MEASl0rdasListl0rdsCoilElementIDdtCoilID,...
          '"','');;
    end;
    if isempty(SeriesInfo(s).SliceThickness) & ...
        isfield(mdhStruct,'sSliceArraydasSlicel1rddThickness')
      SeriesInfo(s).SliceThickness = ...
        mdhStruct.sSliceArraydasSlicel1rddThickness;
    end;
    if isfield(mdhStruct,'sSliceArraydasSlicel1rddPhaseFOV')
      SeriesInfo(s).phaseFOV = mdhStruct.sSliceArraydasSlicel1rddPhaseFOV;
    end;
    if isfield(mdhStruct,'sSliceArraydasSlicel1rddReadoutFOV')
      SeriesInfo(s).readoutFOV = mdhStruct.sSliceArraydasSlicel1rddReadoutFOV;
    end;
    if isfield(mdhStruct,'sKSpacedlPhaseEncodingLines')
      SeriesInfo(s).NumberOfPhaseEncodingSteps = mdhStruct.sKSpacedlPhaseEncodingLines;
    end;
    if isfield(mdhStruct,'sKSpacedlBaseResolution')
      SeriesInfo(s).nc = mdhStruct.sKSpacedlBaseResolution;
    end;
    if isfield(mdhStruct,'sKSpacedlPhaseEncodingLines')
      SeriesInfo(s).nr = mdhStruct.sKSpacedlPhaseEncodingLines;
    end;
    if isfield(mdhStruct,'sAdjDatadsAdjVolumeddInPlaneRot')
      SeriesInfo(s).InPlaneRot = mdhStruct.sAdjDatadsAdjVolumeddInPlaneRot;
    elseif isfield(mdhStruct,'sSliceArraydasSlicel1rddInPlaneRot')
      SeriesInfo(s).InPlaneRot = mdhStruct.sSliceArraydasSlicel1rddInPlaneRot;
    end;
    if isfield(mdhStruct,'sSliceArraydlSize')
      SeriesInfo(s).nslices = mdhStruct.sSliceArraydlSize;
    end;
    if isfield(mdhStruct,'sDiffusiondlDiffDirections')
      SeriesInfo(s).ndiffdirs = mdhStruct.sDiffusiondlDiffDirections;
    end;
    if isfield(mdhStruct,'sDiffusiondalBValuel1r')
      SeriesInfo(s).bval = mdhStruct.sDiffusiondalBValuel1r;
    end;
    if ~isfield(SeriesInfo(s).info, 'Private_0019_100c')
      if isfield(mdhStruct,'sWiPMemBlockdalFreel8r')
        SeriesInfo(s).nb0 = mdhStruct.sWiPMemBlockdalFreel8r;
      elseif isfield(mdhStruct,'sWiPMemBlockdalFreel1r') % is this correct?
        SeriesInfo(s).nb0 = mdhStruct.sWiPMemBlockdalFreel1r;
      end;
      if isfield(mdhStruct,'sDiffusiondalBValuel1r')
        SeriesInfo(s).bval = mdhStruct.sDiffusiondalBValuel1r;
      end;
    end;
    if strcmp(SeriesInfo(s).MRAcquisitionType,'3D')
      pos = strfind(txt,'sKSpace.lImagesPerSlab');
      n = regexp(txt(pos(1):end),...
                 'sKSpace.lImagesPerSlab\s*=\s*(?<slices>\d+)','names');
      SeriesInfo(s).ImagesInAcquisition = str2num(n(1).slices);
    else
      pos1 = strfind(txt,'sSliceArray.lSize');
      n1 = regexp(txt(pos1(1):end),...
                  'sSliceArray.lSize\s*=\s*(?<slices>\d+)','names');
      pos2 = strfind(txt,'lContrasts');
      n2 = regexp(txt(pos2(1):end),'lContrasts\s*=\s*(?<contrasts>\d+)','names');
      if length(n1)>0 & length(n2)>0
        SeriesInfo(s).ImagesInAcquisition = ...
          str2num(n1(1).slices) * str2num(n2(1).contrasts);
      end
    end
  else
    fprintf('%s: WARNING: DICOM tag Private_0029_1020 (mdh) not found\n',...
      mfilename);
  end
  SeriesInfo(s).magnitude_flag = ...
    ~isempty(findstr('ORIGINAL\PRIMARY\M\ND',SeriesInfo(s).ImageType)) | ...
    ~isempty(findstr('ORIGINAL\PRIMARY\OTHER',SeriesInfo(s).ImageType));
  % set flags to determine phase-encode direction for EPI scans
  pe_tag_valid = false;   % phase encoding tag
  if ~isempty(cell2mat(regexp(SeriesInfo(s).SequenceName,...
                             {'ep_b','epse2d1','epfid2d1'})))
    if ~isempty(SeriesInfo(s).InPlaneRot)
      pe_tag_valid = logical(1);
       if abs(SeriesInfo(s).InPlaneRot)>2.0
          SeriesInfo(s).pe_rev_flag = pe_tag_valid;
          SeriesInfo(s).pe_for_flag = logical(0);
       else
          SeriesInfo(s).pe_for_flag = pe_tag_valid;
          SeriesInfo(s).pe_rev_flag = logical(0);
       end;
    else
      pe_tag_valid = isfield(SeriesInfo(s).info,'Private_0051_100b') &&...
                   ischar(SeriesInfo(s).info.Private_0051_100b);
      SeriesInfo(s).pe_rev_flag = pe_tag_valid &&...
                  ~isempty(regexp(SeriesInfo(s).info.Private_0051_100b,'s'));
      SeriesInfo(s).pe_for_flag = pe_tag_valid &&...
                  isempty(regexp(SeriesInfo(s).info.Private_0051_100b,'s'));
    end
  end;
  SeriesInfo(s).pe_extra_val = mmil_getfield(SeriesInfo(s),'Private_0019_1016');
end;

