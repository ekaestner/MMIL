function B1corr_flag = mmil_check_B1corr(ScanInfo,scantypes,scannum)
%function B1corr_flag =mmil_check_B1corr(ScanInfo,scantypes,scannum)
%
% Purpose: Checks if B1 correction is needed 
%         (for at least one scan in ScanInfo)
%
% ScanInfo: struct array containing info about each T1-weighted scan
%     from ContainerInfo
%
% scantypes: cell array of scan types to correct
%
% scannum : scan number in the list of scantypes
% 
% Created:  04/12/12 by Vijay Venkatraman
% Last Mod: 02/28/14 by Don Hagler
%

if ~mmil_check_nargs(nargin,2), return; end;

if ~iscell(scantypes), scantypes = {scantypes}; end;
if ~exist('scannum','var'), scannum = []; end;
B1corr_flag = 0;
for f=1:length(scantypes)
  fnamestem = scantypes{f};
  tmpinfo_array = mmil_getfield(ScanInfo,fnamestem,[]); 
  if ~isempty(scannum)
    tmpinfo_array = tmpinfo_array(scannum);
  end;
  for i=1:length(tmpinfo_array)
    tmpinfo = tmpinfo_array(i);
    PUREflag = mmil_getfield(tmpinfo,'PUREflag',1);
    coilname = mmil_getfield(tmpinfo,'ReceiveCoilName','HEAD');
    coilflag = ismember(coilname,{'HEAD','BODY'});
    if ~coilflag && ~PUREflag
      % B1 correction needed if any scan is not head coil or any scan has PUREflag=0
      B1corr_flag = 1;
      break;
    end;
  end;
end;

return;

