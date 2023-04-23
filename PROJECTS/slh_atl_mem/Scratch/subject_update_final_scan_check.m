clear; clc;

out_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/slh_atl_mem/Analysis/subject_check';

scn_dir = '/home/ekaestne/PROJECTS/OUTPUT/slh_atl_mem/Data';

%% Load
% Subjects
fcfg = [];
fcfg.dta_loc = [ out_dir '/' 'Rebecca_final_demographics.csv'];
[ sbj_dta, sbj_sbj, sbj_col ] = ejk_dta_frm(fcfg);

% QC
fcfg = [];
fcfg.dta_loc = [ out_dir '/' 'QC_overall_ejk.csv'];
[ qal_dta, qal_sbj, qal_col ] = ejk_dta_frm(fcfg);

% DTI
fcfg = [];
fcfg.dta_loc = [ scn_dir '/' 'FinalSample.csv'];
[ dti_dta, dti_sbj, dti_col ] = ejk_dta_frm(fcfg);

%% Go through
sbj_scn_sbj = sbj_sbj;
sbj_scn_col = { 'has_mri' 'has_dti' 'class' 'surgery' 'rebecca_notes' 'rebecca_review' 'dan_meeting_notes' 'dan_email_notes' 'misc_notes' };
sbj_scn_dta = cell( numel(sbj_sbj), numel(sbj_scn_col) );
for iS = 1:numel(sbj_scn_sbj)
    
    if ~sum(strcmpi(dti_sbj,sbj_scn_sbj{iS}))==0
        sbj_scn_dta{iS,1} = ~isempty(dti_dta{strcmpi(dti_sbj,sbj_scn_sbj{iS}),strcmpi(dti_col,'xLeft_Hippocampus')});
        sbj_scn_dta{iS,2} = ~isempty(dti_dta{strcmpi(dti_sbj,sbj_scn_sbj{iS}),strcmpi(dti_col,'xL_Unc')});
    else
        if ~sum(strcmpi(qal_sbj,sbj_scn_sbj{iS}))==0
            error('line 34')
        else
            sbj_scn_dta{iS,1} = 0;
            sbj_scn_dta{iS,2} = 0;
        end
    end
    sbj_scn_dta{iS,3} = sbj_dta{iS, strcmpi(sbj_col,'Class')};
    sbj_scn_dta{iS,4} = sbj_dta{iS, strcmpi(sbj_col,'SurgeryCategory')};
    sbj_scn_dta{iS,5} = sbj_dta{iS, strcmpi(sbj_col,'RebeccaSurgery')};
    sbj_scn_dta{iS,6} = sbj_dta{iS, strcmpi(sbj_col,'RebeccaReview')};
    sbj_scn_dta{iS,7} = sbj_dta{iS, strcmpi(sbj_col,'DanMeetingNotes')};
    sbj_scn_dta{iS,8} = sbj_dta{iS, strcmpi(sbj_col,'DanEmailNotes')};
    sbj_scn_dta{iS,9} = sbj_dta{iS, strcmpi(sbj_col,'Notes')};
    
end

%% Save out
cell2csv([ scn_dir '/' 'Rebecca_Scan_Check.csv'], [ 'sbj' sbj_scn_col ; sbj_scn_sbj sbj_scn_dta ]);