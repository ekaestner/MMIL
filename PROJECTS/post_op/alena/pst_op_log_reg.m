clear; clc;

aln_dta = mmil_readtext('/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena/LeftTLEs_LM_forAUC.csv');
    reo_dta = cell2mat(aln_dta(2:end,2));
    reo_dta_cll(reo_dta==0) = {'NoDecline'};
    reo_dta_cll(reo_dta==1) = {'Decline'};
    aln_dta(2:end,2) = reo_dta_cll';
    
% L-TLE
lcfg = [];

lcfg.sbj_grp = aln_dta( :, 1:2);
lcfg.lbl_ord = { 'Decline' 'NoDecline' };

lcfg.dta     = aln_dta( 2:end, [1 3:6]);
lcfg.dta_lbl = aln_dta( 1, [1 3:6]);

lcfg.mdl     = { { 'PRE_Logical_Memory_II_SS' 'hippocampal_volume_LI' } ...
                 ...{ 'FA_Unc_LI' 'FA_Entorhinal_LI' } ...
                 { 'PRE_Logical_Memory_II_SS' 'hippocampal_volume_LI' 'FA_Unc_LI' 'FA_Entorhinal_LI' } };
lcfg.mdl_nme = { 'Clinical' ...
                 ...'WhiteMatter' ...
                 'Clinical+WhiteMatter' };

lcfg.nrm_grp = 0;

lcfg.mdl_cmp_plt = { [1 2] };
lcfg.mld_cmp_col = { { rgb('green') rgb('purple') } };
lcfg.mdl_cmp_nme = { [ 'LMII' '_' 'Logistic'] };

lcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena/LogMemII/';

mmil_log_reg(lcfg)










