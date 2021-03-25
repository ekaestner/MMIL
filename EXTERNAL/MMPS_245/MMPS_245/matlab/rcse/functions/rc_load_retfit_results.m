function retfit_results = rc_load_retfit_results(retfit_stem)
%function retfit_results = rc_load_retfit_results(retfit_stem)
%
% Required Parameters:
%   retfit_stem: full path retfit file stem
%     e.g. 'retfit_path/retfit'
%
% Created:  02/06/09 by Don Hagler
% Last Mod: 10/25/13 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

tags = {'data','fit_data','area_masks','area_data'};

retfit_results = [];
hemilist = {'lh','rh'};
for h=1:length(hemilist)
  hemi = hemilist{h};
  fname = sprintf('%s_results-%s.mat',retfit_stem,hemi);
  if ~exist(fname,'file'), error('file %s not found',fname); end;
  results = load(fname);
  for i=1:length(tags)
    tag = tags{i};
    retfit_results.([hemi '_' tag]) = mmil_getfield(results,tag);
  end;
end;

