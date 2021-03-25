clear; clc

%% Amplitude - VISUAL VS AUDITORY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HEMISPHERE
% onset
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'hem';
fcfg.loc_nme = {  'vis_act' 'aud_act' 'bim_act' };

fcfg.col      = [ 1         2        3 ];
fcfg.ylm_typ  = { 'maxmin'  'maxmin' 'maxmin' };

fcfg.hms     = {'lhs' 'rhs'};
fcfg.plt     = 1;
fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure6' '/' 'Table_amplitude' '/'];
hms_sig_eff  = mmil_ranksum_region_test(fcfg);

hms_sig_eff.cmb_reg

% WITHIN REGIONS
% onset
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'wth_reg';
fcfg.loc_nme = { 'vis_act' 'aud_act' 'bim_act'};
fcfg.use_col = { [1 2]     };

fcfg.hms     = {'lhs' 'lhs'};
fcfg.plt     = 1;
fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure6' '/' 'Table_amplitude' '/'];
wth_sig_lhs_eff = mmil_ranksum_region_test(fcfg);

wth_sig_lhs_eff.cmb_reg.pvl

% BETWEEN Regions
% amplitude
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = { 'vis_act' 'aud_act' 'bim_act' };

fcfg.col      = [ 1          2         3 ];
fcfg.ylm_typ  = { 'maxmin'   'maxmin'  'maxmin' };

fcfg.hms        = {'lhs' 'lhs'};
fcfg.plt        = 1;
fcfg.out        = [fcfg.clr_fld '/' 'manuscript' '/' 'figure6' '/' 'Table_amplitude' '/'];
bet_sig_lhs_eff = mmil_ranksum_region_test(fcfg);


% spc_tst = [ bet_sig_lhs_eff.cmb_reg.vis_act.pvl{9,1} bet_sig_lhs_eff.cmb_reg.vis_act.pvl{1,6}  bet_sig_lhs_eff.cmb_reg.vis_act.pvl(9,6) ; ...
%             bet_sig_lhs_eff.cmb_reg.vis_act.pvl{9,1} bet_sig_lhs_eff.cmb_reg.vis_act.pvl{1,8}  bet_sig_lhs_eff.cmb_reg.vis_act.pvl(9,8) ; ...
%             bet_sig_lhs_eff.cmb_reg.vis_act.pvl{9,1} bet_sig_lhs_eff.cmb_reg.vis_act.pvl{1,10} bet_sig_lhs_eff.cmb_reg.vis_act.pvl(9,10) ]
%           [bet_sig_lhs_eff.cmb_reg.vis_act.pvl([1 end],:) ]
% 
% spc_tst = [ bet_sig_lhs_eff.cmb_reg.aud_act.pvl{9,1} bet_sig_lhs_eff.cmb_reg.aud_act.pvl{1,6}  bet_sig_lhs_eff.cmb_reg.aud_act.pvl(9,6) ; ...
%             bet_sig_lhs_eff.cmb_reg.aud_act.pvl{9,1} bet_sig_lhs_eff.cmb_reg.aud_act.pvl{1,8}  bet_sig_lhs_eff.cmb_reg.aud_act.pvl(9,8) ; ...
%             bet_sig_lhs_eff.cmb_reg.aud_act.pvl{9,1} bet_sig_lhs_eff.cmb_reg.aud_act.pvl{1,10} bet_sig_lhs_eff.cmb_reg.aud_act.pvl(9,10) ]
%           [bet_sig_lhs_eff.cmb_reg.aud_act.pvl([1 end],:) ]

spc_tst = [ bet_sig_lhs_eff.cmb_reg.bim_act.pvl{9,1} bet_sig_lhs_eff.cmb_reg.bim_act.pvl{1,6}  bet_sig_lhs_eff.cmb_reg.bim_act.pvl(9,6) ; ...
            bet_sig_lhs_eff.cmb_reg.bim_act.pvl{9,1} bet_sig_lhs_eff.cmb_reg.bim_act.pvl{1,8}  bet_sig_lhs_eff.cmb_reg.bim_act.pvl(9,8) ; ...
            bet_sig_lhs_eff.cmb_reg.bim_act.pvl{9,1} bet_sig_lhs_eff.cmb_reg.bim_act.pvl{1,10} bet_sig_lhs_eff.cmb_reg.bim_act.pvl(9,10) ]
          [bet_sig_lhs_eff.cmb_reg.bim_act.pvl([1 end],:) ]

%% Amplitude - Bi-modal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% Bimodal Versus Uni-Modal
% amplitude
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'wth_reg';
fcfg.loc_nme = { '' '' '' 'vis_uni' 'vis_bim' 'tot_uni' 'aud_uni' 'aud_bim' };
fcfg.use_col = { [4 5] [6 5] [7 8] };

fcfg.hms     = {'lhs' 'lhs'};
fcfg.plt     = 1;
fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure6' '/' 'Table_amplitude' '/'];
wth_sig_lhs_eff = mmil_ranksum_region_test(fcfg);

wth_sig_lhs_eff.cmb_reg.pvl
% wth_sig_lhs_eff.tot_reg.pvl