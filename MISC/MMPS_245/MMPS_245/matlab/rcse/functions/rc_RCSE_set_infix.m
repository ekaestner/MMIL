function infix = rc_RCSE_set_infix(parms)
%function infix = rc_RCSE_set_infix(parms)
%
% Purpose: create string to indicate values of some parameters
%
% Created:  03/05/11 by Don Hagler
% Last Mod: 03/05/11 by Don Hagler
%
%

if ~mmil_check_nargs(nargin,1), return; end;

infix = [];

% channel types
if parms.usegrad_flag, infix = [infix 'grad']; end;
if parms.usemag_flag, infix = [infix 'mag']; end;
if parms.useEEG_flag, infix = [infix 'EEG']; end;

% forward model
if parms.bem_flag
  infix = [infix '_bem'];
else
  infix = [infix '_sph'];
end;

% constraints
if parms.loose_flag
  infix = sprintf('%s_loose_smf%0.5f_wtang%0.3f',...
    infix,parms.indy_smfact,parms.loose_tang_weight);
elseif parms.indy_locs_flag
  infix = sprintf('%s_indy_smf%0.5f',...
    parms.indy_smfact);
else
  infix = [infix '_full'];
end;
if parms.indy_locs_flag
  infix = sprintf('%s_eccsmf%0.5f',infix,parms.ecc_smfact);
  infix = sprintf('%s_ulsmf%0.1f',infix,parms.upperlower_smfact);
  infix = sprintf('%s_hemismf%0.1f',infix,parms.hemi_smfact);
end;

% extra retinotopic and non-retinotopic dipoles
if parms.ret_dips_weight | parms.nonret_dips_weight
  if parms.ret_dip_qfield_flag & parms.ret_dips_weight~=0
    infix = sprintf('%s_retd%0.2f_qfield_nretd%0.2f',...
      infix,parms.ret_dips_weight,parms.nonret_dips_weight);
  else
    infix = sprintf('%s_retd%0.2f_nretd%0.2f',...
      infix,parms.ret_dips_weight,parms.nonret_dips_weight);
  end;
end;

% SNR
if parms.SNR<1
  infix = sprintf('%s_SNR%0.1f',infix,parms.SNR);
else
  infix = sprintf('%s_SNR%d',infix,parms.SNR);
end;

% noise covariance
if parms.ncov_type==1
  infix = [infix '_avgncov'];
elseif parms.ncov_type==2
  infix = [infix '_rawncov'];
end;

% which conditions
if length(parms.hemivec)==1
  infix = sprintf('%s_hemi%d',infix,parms.hemivec);
end;
if length(parms.uplowvec)==1
  infix = sprintf('%s_uplow%d',infix,parms.uplowvec);
end;
if length(parms.eccvec)==1
  infix = sprintf('%s_ecc%d',infix,parms.eccvec);
elseif ~isempty(parms.eccvec)
  infix = sprintf('%s_ecc%s',...
    infix,sprintf('_%d',parms.eccvec));
end;
if length(parms.contvec)==1
  infix = sprintf('%s_cont%d',...
    infix,parms.contvec);
elseif ~isempty(parms.contvec)
  infix = sprintf('%s_cont%s',...
    infix,sprintf('_%d',parms.contvec));
end;

% which visual areas
if ~isempty(parms.use_areas)
  infix = sprintf('%s_areas%s',...
    infix,sprintf('_%d',parms.use_areas));
end;

% input cond_info file with dipole offsets
if parms.cond_offsets_flag
  infix = sprintf('%s_condoffsets',infix);
end;

% restrict dipole patches
if parms.restrict_hemi_flag
  infix = sprintf('%s_rxhemi',infix);
end;
if parms.restrict_uplow_flag
  infix = sprintf('%s_rxuplow',infix);
end;

% for offset search
if parms.polarity_penalty~=0
  infix = sprintf('%s_polpenalty%0.0f',...
    infix,parms.polarity_penalty);
end;
if parms.fit_range_flag~=0
  infix = [infix '_fitrange'];
end;
if parms.cond_offsets_flag==2
  infix = [infix '_EEG'];
end;
if parms.offset_niters>0
  infix = sprintf('%s_offsetgroupconds%dpatches%dareas%d',...
    infix,parms.offset_group_flag,...
    parms.offset_group_patches_flag,parms.offset_group_areas_flag);
end;
if ~isempty(parms.offset_const_areas)
  infix = sprintf('%s_constoffareas%s',...
    infix,sprintf('_%d',parms.offset_const_areas));
end;
