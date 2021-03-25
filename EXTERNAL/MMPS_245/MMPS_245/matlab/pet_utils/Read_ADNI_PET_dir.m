function [N, pibflag] = Read_ADNI_PET_dir(filedir)
%function [N, pibflag] = Read_ADNI_PET_dir(filedir)
%
% Last Mod: 04/06/10 by Don Hagler
%

N = [];
pibflag = 0;

% Look for ECAT format
dlist = dir(sprintf('%s/*.v',filedir));
if length(dlist)>0
  if length(dlist)>1
    fprintf(1,'%s: More than one *.v files in dir %s\n',mfilename,filedir);
    return
  end
  fname = sprintf('%s/%s',filedir,dlist(1).name);
  N = ReadECATfile(fname);
  if ~isempty(regexp(N{1}.extras.mh.ISOTOPE_NAME,'C-11')) || ...
      ~isempty(regexp(N{1}.extras.mh.RADIOPHARAMCEUTICAL,'C-11'))
    pibflag = 1;
  end
  return
end

% Look for Interfile format
dlist = dir(sprintf('%s/*.i.hdr',filedir));
if length(dlist)>0
  for i = 1:length(dlist)
    fnames{i} = sprintf('%s/%s',filedir,dlist(i).name);
  end
  N = ReadInterfiles(fnames);
  if ~isempty(regexp(N{1}.info.DoseType,'C-11'))
    pibflag = 1;
  end
  return
end

% Look for DICOM format
dlist = dir(sprintf('%s/*.dcm',filedir));
if length(dlist)>0
  for i = 1:length(dlist)
    fnames{i} = sprintf('%s/%s',filedir,dlist(i).name);
  end
  N = ReadPETdcmfiles(fnames);
  info = dicominfo(fnames{1});
  try
    if isfield(info.RadiopharmaceuticalInformationSequence.Item_1,'Radiopharmaceutical') && ...
        ~isempty(regexp(info.RadiopharmaceuticalInformationSequence.Item_1.Radiopharmaceutical,'PIB'))
      pibflag = 1;
    elseif isfield(info.RadiopharmaceuticalInformationSequence.Item_1.RadionuclideCodeSequence.Item_1,'CodeMeaning') && ...
        (~isempty(regexp(info.RadiopharmaceuticalInformationSequence.Item_1.RadionuclideCodeSequence.Item_1.CodeMeaning,'C.11')) || ...
        ~isempty(regexp(info.RadiopharmaceuticalInformationSequence.Item_1.RadionuclideCodeSequence.Item_1.CodeMeaning,'PIB')))
      pibflag = 1;
    end
  catch
    fprintf('%s: WARNING: missing fields in dicom header:\n%s\n',...
      mfilename,lasterr);
  end;
  return
end

fprintf(1,'%s: No PET files found in dir %s\n',mfilename,filedir);
