function [mdhStruct,txt,errcode] = mmil_read_mdhStruct(dcminfo)
%function [mdhStruct,txt,errcode] = mmil_read_mdhStruct(dcminfo)
%
% Created:  08/31/11 by Don Hagler
% Prev Mod: 03/31/15 by Don Hagler
% Prev Mod: 01/03/17 by Don Hagler
% Last Mod: 05/02/17 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
mdhStruct = []; txt = []; errcode = 0;

tag = 'Private_0029_1020';
startstr = '### ASCCONV BEGIN';
markstr = '###';
endstr = '### ASCCONV END ###';

if isfield(dcminfo,tag)
  txt = char(dcminfo.(tag)');
  try
    % find starting position
    startpos = strfind(txt,startstr);
    startpos = startpos(1) + length(startstr);
    % jump to end of starting comment
    markpos = strfind(txt(startpos:end),markstr);
    startpos = startpos + markpos(1) + length(markstr);
    % find ending position
    endpos = strfind(txt,endstr)-1;
    % exclude end points prior to start points
    endpos = endpos(endpos>startpos);
    % select string between start point and first end point    
    hdrstr = txt(startpos:endpos(1));
    warning off;
    mdhStruct = read_mdhStruct_from_string(hdrstr);
    warning on;
  catch
    fprintf('%s: WARNING: failed to read mdhStruct from dicom header string\n',...
      mfilename);
    errcode = 1;
  end;
end;

