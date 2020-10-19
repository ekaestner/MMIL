function wforms = rc_RCSE_extract_source_wforms(parms,avg_data,results)
%function wforms = rc_RCSE_extract_source_wforms(parms,avg_data,results)
%
% Purpose: extract source waveforms from results struct
%
% Required Input:
%   parms: RCSE parms struct
%   avg_data: average data struct
%   results: results struct
%
% Created:  03/04/11 by Don Hagler
% Last Mod: 04/27/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,3), return; end;

fprintf('%s: extracting source waveforms...\n',mfilename);
retmap = parms.best_retmap;
num_areas = retmap.num_areas;
if ~parms.indy_locs_flag
  num_sources = 1;
else
  num_sources = retmap.num_locs;
end;  
num_tpoints = length(avg_data.averages(1).time);
num_locs = retmap.num_locs;
wforms = zeros(num_areas,num_sources,num_tpoints,parms.ncontrasts);
for c=1:parms.ncontrasts
  for a=1:num_areas
    if parms.loose_flag
      j = 1+(a-1)*num_locs*3;
      k = j+num_locs*3-1;
      srange = [j:3:k];
    elseif parms.indy_locs_flag
      j = 1+(a-1)*num_locs;
      k = j+num_locs-1;
      srange = [j:k];
    else
      srange = a;
    end
    wforms(a,:,:,c) = (results.S(:,srange,c))';
  end
end;

