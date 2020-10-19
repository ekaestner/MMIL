function REC_MMIL_Status_Report(ProjID)
%function REC_MMIL_Status_Report(ProjID)
%
% Created:   04/14/09 Alain Koyama
% Last Mod:  11/08/12 Don Hagler
%

if ~exist('ProjID','var'), ProjID = []; end;

if isempty(ProjID)
   error('Empty ProjID');
end

[RootDirs,ProjInfo] = REC_MMIL_RootDirs(ProjID);

analysis_outdir = 'analysis';

for p=1:length(ProjInfo)
  [RootDirs,ProjInfo] = REC_MMIL_RootDirs(ProjID);
  if isempty(RootDirs.orig_pet)
    pet_flag = 0;
  else
    pet_flag = 1;
  end
  MMIL_Status_Report(RootDirs,...
    'ProjID',ProjID,...
    'analysis_outdir',analysis_outdir,...
    'pet_flag',pet_flag)
end
