function weights = rc_create_gaussian_weights(surf,pol_r,pol_i,ecc_r,ecc_i,...
  areamask,hemi,stim_pol_phases,stim_ecc_phases,stim_pol_size,stim_ecc_size,...
  pol_sigma,ecc_sigma,logtrans_flag,ecc_maxvang,ecc_minvang);
%function weights = rc_create_gaussian_weights(surf,pol_r,pol_i,ecc_r,ecc_i,...
%  areamask,hemi,stim_pol_phases,stim_ecc_phases,stim_pol_size,stim_ecc_size,...
%  pol_sigma,ecc_sigma,[logtrans_flag],[ecc_maxvang],[ecc_minvang]);
%
% Required Input:
%   surf: freesurfer surface structure containing:
%     nverts: number of vertices
%     nfaces: number of faces (triangles)
%     faces:  vertex numbers for each face (3 corners)
%     coords: x,y,z coords for each vertex
%     (see fs_load_subj and fs_read_surf)
%  pol_r: real component of polar angle map
%  pol_i: imaginary component of polar angle map
%  ecc_r: real component of eccentricity map
%  ecc_i: imaginary component of eccentricity map
%    these should be vectors of values for each vertex (can be sparse)
%    (see fs_read_wfile)
%  areamask: vector of vertex numbers to mask a single visual area
%  hemi: cortical hemisphere ('lh' or 'rh')
%  stim_pol_phases: vector of phases corresponding to visual field polar angles
%  stim_ecc_phases: vector of phases corresponding to visual field eccentricities
%    note: number of stim_pol_phases must match number of stim_ecc_phases
%    these phases should be in units of cycles (0->1)
%  stim_pol_size: polar angle size of stimulus (in cycles 0->1)
%  stim_ecc_size: eccentricity size of stimulus (in cycles 0->1)
%  pol_sigma: polar angle receptive field sigma (gaussian blur)
%  ecc_sigma: eccentricity receptive field sigma (gaussian blur)
%  logtrans_flag: [0|1] assume log transform relation between phase and eccentricity
%    {default = 1}
%  ecc_minvang: minimum eccentricity visual angle (deg) (used for logtrans)
%    {default = 0.2}
%  ecc_maxvang: maximum eccentricity visual angle (deg) (used for logtrans)
%    {default = 12.5}
%
% Output:
%   weights: sparse matrix containing weights for each pol/ecc phase and each vertex
%
% Created: 12/05/06 Don Hagler
% Lst Mod: 02/19/11 Don Hagler
%

%% todo: use varargin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,10), return; end;
if ~exist('logtrans_flag','var'), logtrans_flag = 1; end;
if ~exist('ecc_minvang','var'), ecc_minvang = 0.2; end;
if ~exist('ecc_maxvang','var'), ecc_maxvang = 12.5; end;
weights = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check for errors
if ~ismember(hemi,{'lh','rh'})
  fprintf('%s: error: hemi must be ''lh'' or ''rh''\n',mfilename);
  return;
end;

if length(stim_pol_phases) ~= length(stim_ecc_phases)
  fprintf('%s: error: number of stim_pol_phases must match number of stim_ecc_phases\n',...
    mfilename);
  return;
end;
if length(find(stim_pol_phases>1))>0 | length(find(stim_pol_phases<-1))>0
  fprintf('%s: error: stim_pol_phases must be between -1 and 1\n',...
    mfilename);
  return;
end;
if length(find(stim_ecc_phases>1))>0 | length(find(stim_ecc_phases<-1))>0
  fprintf('%s: error: stim_ecc_phases must be between -1 and 1\n',...
    mfilename);
  return;
end;
if stim_pol_size>1 | stim_pol_size<0
  fprintf('%s: error: stim_pol_size must be between 0 and 1\n',...
    mfilename);
  return;
end;
if stim_ecc_size>1 | stim_ecc_size<0
  fprintf('%s: error: stim_ecc_size must be between 0 and 1\n',...
    mfilename);
  return;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% apply mask to data
pol_r = full(pol_r(areamask));
pol_i = full(pol_i(areamask));
ecc_r = full(ecc_r(areamask));
ecc_i = full(ecc_i(areamask));

% calculate phase in units of cycles
pol_phase = atan2(pol_i,pol_r)/(2*pi);
ecc_phase = atan2(ecc_i,ecc_r)/(2*pi);

% add 1 to negative ecc_phase
tmp = find(ecc_phase<0);
ecc_phase(tmp)=ecc_phase(tmp)+1;

nstims = length(stim_pol_phases);
nverts = surf.nverts;
weights = sparse(nstims,nverts);

if stim_pol_size>0 | stim_ecc_size>0
  for i=1:nstims
    th1 = stim_pol_phases(i) - stim_pol_size;
    th2 = stim_pol_phases(i) + stim_pol_size;
    r1 =  stim_ecc_phases(i) - stim_ecc_size;
    r2 =  stim_ecc_phases(i) + stim_ecc_size;
    if logtrans_flag
      % log transform of r1 and r2
      r1 = log(r1*ecc_maxvang/ecc_minvang)/log(ecc_maxvang/ecc_minvang);
      r2 = log(r2*ecc_maxvang/ecc_minvang)/log(ecc_maxvang/ecc_minvang);
    end;
    a = 0.5*pol_sigma*ecc_sigma*pi;
    b = sqrt(2)/(2*pol_sigma);
    c = sqrt(2)/(2*ecc_sigma);
    % this function is the result of integrating over the range of r and theta
    weights(i,areamask) = a*(erf(b*(pol_phase-th2))-erf(b*(pol_phase-th1))).*...
                            (erf(b*(ecc_phase-r2))-erf(b*(ecc_phase-r1)));
    % repeat for pol_phases + 1
    weights(i,areamask) = weights(i,areamask)' +...
                          a*(erf(b*(pol_phase-th2+1))-erf(b*(pol_phase-th1+1))).*...
                          (erf(b*(ecc_phase-r2))-erf(b*(ecc_phase-r1)));
  end;
  % normalze weights
  maxw = max(weights,[],2);
  weights = weights ./ (maxw*ones(1,size(weights,2))+eps);
else
  for i=1:nstims
    th = stim_pol_phases(i);
    r = stim_ecc_phases(i);
    weights(i,areamask) = exp(-0.5*(((pol_phase-th)/pol_sigma).^2 +...
                                    ((ecc_phase-r)/ecc_sigma).^2));
    % repeat for pol_phases + 1
    weights(i,areamask) = weights(i,areamask)' +...
                          exp(-0.5*(((pol_phase-th+1)/pol_sigma).^2 +...
                                    ((ecc_phase-r)/ecc_sigma).^2));
  end;
end;


return;

