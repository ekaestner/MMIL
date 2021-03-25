function [arg_groups,extra_tags,strip_flag] = mmil_args_rc_fit_wforms()
%function [arg_groups,extra_tags,strip_flag] = mmil_args_rc_fit_wforms()
%
% Purpose: for use with MMIL_Args in selecting context relevant parameters
%   from parameter struct arrays (e.g. parms, ProjInfo)
%
% Output:
%   arg_groups: cell array of fieldname prefixes (e.g. 'BOLD', 'PROC', etc.)
%     specifying group of parameter names
%   extra_tags: cell array of parameter names
%   strip_flag: [0|1] whether arg_group prefixes should be stripped from
%     final parameter names
%
% Created:  05/13/11 by Don Hagler
% Last Mod: 07/22/11 by Don Hagler
%

arg_groups = [];
extra_tags = {'condition_valuse' 'area_names' 'area_colors' 'niters' 'search_type'...
  'stepsize' 'rand_init_flag' 'ncomponents' 'polarity' 'latency' 'amplitude'...
  'rise_tc' 'fall_tc' 'latency_bounds' 'amplitude_bounds' 'rise_tc_bounds'...
  'fall_tc_bounds' 'sfreq' 't0' 't1' 'time' 'contrast_latency_flag' 'outdir'...
  'outstem' 'delay_sf' 'DiffMinChange' 'TolFun' 'TolX' 'plotflag' 'xlim'...
  'ylim' 'use_areas' 'linewidth_data' 'linewidth_fit' 'visible_flag'...
};
strip_flag = 0;

