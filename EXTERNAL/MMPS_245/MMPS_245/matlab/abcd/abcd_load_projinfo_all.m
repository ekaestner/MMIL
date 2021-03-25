function projinfo = abcd_load_projinfo_all(parms)
%function abcd_load_projinfo_all(parms)
%
% Optional input:
%
% Created:  05/18/17 by Feng Xue
% Last Mod: 05/18/17 by Feng Xue
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if ~exist(parms.fname_projinfo,'file'), error('info file %s not found',parms.fname_projinfo); end;
  tmp = mmil_csv2struct (parms.fname_projinfo);
  tmp = mmil_sortstruct (tmp,{'ProjID','sourceID'});
  idx = ~cellfun('isempty',strfind({tmp.ProjID},parms.instem));
  projinfo = tmp(idx);
return;
