function [Channel,CoilType] = abcd_get_coilinfo(dcminfo,vendor)
%function [Channel,CoilType] = abcd_get_coilinfo(dcminfo,vendor)
%
% Purpose: Extract chennel and coiltype based on dicominfo and vendor.
%
% Created:  06/22/17 by Feng Xue
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Channel = [];
  CoilType = [];
  switch vendor
    case 'Siemens'
      [mdhStruct,txt] = mmil_read_mdhStruct(dcminfo);
      if isfield(mdhStruct,'sCoilSelectMeasdaRxCoilSelectDatal0rdasListd__attribute__dsize'), Channel = mdhStruct.sCoilSelectMeasdaRxCoilSelectDatal0rdasListd__attribute__dsize; end;
      index = findstr(txt,'sCoilSelectMeas.aRxCoilSelectData[0].asList[0].sCoilElementID.tCoilID');
      if ~isempty(index)
        txt = txt(index:end);
        index = findstr(txt,'"');
        if length(index) >4, CoilType = txt(index(2)+1:index(3)-1); end;
      else
        if isfield(mdhStruct,'sCoilSelectMeasdaRxCoilSelectDatal0rdasListl0rdsCoilElementIDdt'), CoilType = regexprep(mdhStruct.sCoilSelectMeasdaRxCoilSelectDatal0rdasListl1rdsCoilElementIDdt,'"',''); end;
      end
    case 'GE'
      CoilType = dcminfo.ReceiveCoilName;
      Channel = regexprep(CoilType,'[^0-9]*(\d*)[^0-9]*','$1');
    case 'Philips'
      CoilType = dcminfo.ReceiveCoilName;
      Channel = regexprep(CoilType,'[^0-9]*(\d*)[^0-9]*','$1');
  end
return
