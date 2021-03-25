function avg_angle = rc_nbrhd_dispersion(v,surf,dip_info,nsteps)
%function avg_angle = rc_nbrhd_dispersion(v,surf,dip_info,[nsteps])
%
% Purpose: calculate variation of dipole orientation
%   as a function of neighbor distance from vertex v
%   returns the average angle of the normal vector of v and its neighbors
%
% Required Input:
%   v: reference vertex
%   surf: fsurf_tools surf structure
%     (see fs_load_subj)
%   dip_info: dipole info (6 x nverts)
%     (see ts_dip_info and ts_read_dip_file)
% Optional Input:
%   nsteps: neighbor steps from reference vertex
%     {default = 1}
%
% Created:  09/24/07 by Don Hagler
% Last Mod: 01/15/13 by Don Hagler
%

if ~mmil_check_nargs(nargin,3), return; end;
if ~exist('nsteps','var') | isempty(nsteps), nsteps = 1; end;

avg_angle = zeros(nsteps,1);
nvec = squeeze(dip_info(4:6,v));
nvec = nvec / sqrt(sum(nvec.^2));
nbrs = surf.nbrs(v,:);
nbrs = nbrs(nbrs>0);
for n=1:nsteps
  if n>1
    new_nbrs = surf.nbrs(old_nbrs,:);
    new_nbrs = unique(new_nbrs(new_nbrs(:)>0));
    new_nbrs = setdiff(new_nbrs,old_nbrs);
  else
    new_nbrs = nbrs;
  end;
  nvecs = squeeze(dip_info(4:6,new_nbrs));
  avg_dotp = 0;
  for k=1:length(new_nbrs)
    tmp_nvec = nvecs(:,k);
    tmp_nvec = tmp_nvec / sqrt(sum(tmp_nvec.^2));
    dotp = dot(nvec,tmp_nvec);
    avg_dotp = avg_dotp + dotp;
  end;
  avg_dotp = avg_dotp / length(new_nbrs);  
  avg_angle(n) = acos(avg_dotp)*180/pi;
  old_nbrs = new_nbrs;
end;

