function run_dSPM()
%function run_dSPM()
%
% Created:  07/25/11 by Don Hagler
% Last Mod: 08/03/11 by Don Hagler
%

ProjID = 'TEST';
parms = [];
parms.batchname = 'dSPM';
parms.dSPM_prefix = 'dSPM';
parms.dSPM_ncov_type = 2;
parms.dSPM_noisenorm_identity_flag = 0;
parms.dSPMflag = 1;
parms.dSROIflag = 1;
parms.RCSEflag = 0;
parms.dSPM_proc_prefix = 'proc';
parms.dSPM_SNR = 3;
parms.dSPM_usemag_flag = 0;
parms.dSPM_forward_matfile = 'dspm_forward.mat';
parms.dSPM_calc_scalefacts_flag = 0;
parms.dSPM_conditions = [37,73,109];
parms.dSROI_conditions = parms.dSPM_conditions;
parms.dSROI_prefix = parms.dSPM_prefix;
parms.dSROI_outstem = [parms.dSROI_prefix '_ROI_results'];
parms.dSROI_plot_xlim = [-100,300];
parms.dSROI_plot_ylim = [-1,15];
parms.dSROI_legend_flag = 1;
parms.dSROI_roinames = {'DH-pos-V1','DH-dor-occ','DH-lat-occ','DH-ven-occ'};
parms.qcflag = 1;
parms.forceflag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

args = mmil_parms2args(parms);
MMIL_Analyze_MEG_Exams(ProjID,args{:});


