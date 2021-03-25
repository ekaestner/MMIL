function R = rc_calc_ref_absdiff(S,S_ref,varargin)
%function R = rc_calc_ref_absdiff(S,S_ref,[options])
%
% Purpose: calculate overall difference between two waveforms
%
% Required Input:
%   S: source waveform matrix
%   S_ref: reference source waveform matrix
%
% Optional Input:
%  'S_sem_ref': reference source waveform standard error of mean matrix
%     If supplied, will be used to normalize differences between S and S_ref
%     {default = []}
%  'sem_min': set all S_sem_ref values to at least this value
%     {default = 0.5}
%  'srange': vector of start and end samples from which to calculate diff
%    If empty, use all available time samples
%     {default = []}
%  'scale_flag': [0|1] whether to scale S_ref to best match S
%     {default = 0}
%  'indy_wform_flag': [0|1] scale each waveform (areas and contrast) independently
%     Otherwise, scale all based on fit to first area
%     {default = 0}
%
% Created:  03/04/11 by Don Hagler
% Last Mod: 01/19/13 by Don Hagler
%

if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin, { ...
  'S_sem_ref',[],[],...
  'sem_min',0.5,[eps,100],...
  'srange',[],[],...
  'scale_flag',false,[false true],...
  'indy_wform_flag',false,[false true],...
});

if isempty(parms.srange), parms.srange = [1,size(S,1)]; end;

if parms.scale_flag
  fprintf('%s: calculating difference from scaled reference waveforms...\n',mfilename);
  % calculate scaling factors, apply to reference waveforms
  for c=1:size(S,3) % independent scaling factor for each contrast
    if parms.indy_wform_flag % independent scaling factor for each visual area
      for s=1:size(S,2)
        sf = calc_scale_fact(squeeze(S(:,s,c)),squeeze(S_ref(:,s,c)));
        S_ref(:,s,c) = abs(sf) * S_ref(:,s,c);
        if ~isempty(parms.S_sem_ref)
          parms.S_sem_ref(:,s,c) = abs(sf) * parms.S_sem_ref(:,s,c);
        end;
      end;
    else % scaling factor from first visual area applied to all
      sf = calc_scale_fact(squeeze(S(:,1,c)),squeeze(S_ref(:,1,c)));
      S_ref(:,:,c) = abs(sf) * S_ref(:,:,c);
      if ~isempty(parms.S_sem_ref)
        parms.S_sem_ref(:,:,c) = abs(sf) * parms.S_sem_ref(:,:,c);
      end;
    end;
  end;
else
  fprintf('%s: calculating difference from reference waveforms...\n',mfilename);
end;

% calculate difference between waveforms
wform = S(parms.srange(1):parms.srange(2),:);
wform_ref = S_ref(parms.srange(1):parms.srange(2),:);
wform_diff = wform - wform_ref;

% transform difference into z-score
if ~isempty(parms.S_sem_ref)
  wform_sem = parms.S_sem_ref(parms.srange(1):parms.srange(2),:);
  wform_sem = max(parms.sem_min,wform_sem);
  wform_zscore = mean(wform_sem(:)) * wform_diff ./ (wform_sem+eps);
  wform_diff = wform_zscore;
end;

R1 = sum((wform_diff).^2);
R2 = sum((wform.^2 + wform_ref.^2))/2;
R = mean(R1./R2);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sf = calc_scale_fact(v,v_ref)
  sf = v_ref\v;
return;

